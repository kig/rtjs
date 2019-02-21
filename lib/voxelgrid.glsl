// WebGL 2 shader to ray trace.

sampler2D array;

struct Hit {
    float index;
    float closestHit;
};

struct Ray {
    vec3 o;
    vec3 d;
};

vec3 fromGrid(vec3 coord, float scale, vec3 origin) {
    return origin + (coord * scale);
}

vec3 toGrid(vec3 point, float scale, vec3 origin) {
    return (point - origin) / scale;
}

void intersectTri(in Ray ray, in float triIndex, inout Hit closestHit) {
        float off = 4 + triIndex * 18;
        vec3 v0 = readVec3(array, off);
        vec3 e1 = readVec3(array, off+3);
        vec3 e2 = readVec3(array, off+6);

        vec3 h = cross(ray.d, e2);
        float a = dot(e1, h);

        if (a > -0.00000001 && a < 0.00000001) {
            return;
        }

        float f = 1.0 / a;
        vec3 s = sub(ray.o, v0);
        float u = f * dot(s, h);

        if (u < -0.00000001 || u > 1.00000001) {
            return;
        }

        vec3 q = cross(s, e1);
        float v = f * dot(ray.d, q);

        if (v < -0.00000001 || u + v > 1.00000002) {
            return;
        }

        // at this stage we can compute t to find out where
        // the intersection point is on the line
        float t = f * dot(e2, q);

        if (t <= 0.00000001 || t >= closestHit.distance) {
            return;
        }

        closestHit.index = triIndex;
        closestHit.distance = t;
}


void intersectTris(Ray ray, float coff, float childSize, inout Hit closestHit) {
    for (float j = 0; j < childSize; j++) {
        float triIndex = readFloat(array, coff + j);
        if (triIndex < 0) {
            break;
        }
        intersectTri(ray, triIndex, closestHit);
    }
}

/*

class SerializedVG {
    constructor(vgFloat32Array, color) {
        this._color = color;
        this._array = vgFloat32Array;
        this._hit = {
            _color: color,
            _index: 0,
            _e1: vec3(),
            _e2: vec3(),
            _normals: [vec3(), vec3(), vec3()],
            color: function (point) { return this._color },
            normal: function (point) {
                var u = dot(point, this._e1);
                var v = dot(point, this._e2);
                var n0 = this._normals[0];
                var n1 = this._normals[1];
                var n2 = this._normals[2];
                return normalize(mix(mix(n0, n2, v), n1, u));
            }
        };
    }

    intersect(ray) {
        // Step through the voxel grid
        // On encountering a non-negative node,
        // look it up and intersect against it.

        // Do intersect against top-level VG bbox

        // Load _box from _array.
        // Add _box to _array.


        const headOff = this._array[0] * 18 + 4;
        const origin = vec3(
            this._array[headOff],
            this._array[headOff + 1],
            this._array[headOff + 2]
        );
        const size = this._array[headOff + 3];
        const dims = vec3(
            this._array[headOff + 4],
            this._array[headOff + 5],
            this._array[headOff + 6]
        );
        const childCount = this._array[headOff + 7];
        const box = new Box(origin, add(origin, dims));
        const scale = dims.x / size;
        const voxelsOff = headOff + 8;
        const childIndexOff = voxelsOff + size * size * size;
        const childOff = childIndexOff + childCount;

        const t = box.intersect(ray);
        if (t === null) {
            return null;
        }

        // Map hit to voxel coordinates
        const tray = new Ray(ray.o, ray.d, ray.time);
        tray.o = add(tray.o, mulS(tray.d, t.distance + scale * 0.0001));

        const cf = this.toGrid(tray.o, scale, origin);
        const c = floorV(cf);

        const deltaDist = vec3();
        const step = sign(ray.d);
        deltaDist.x = length(mulS(ray.d, 1 / ray.d.x));
        deltaDist.y = length(mulS(ray.d, 1 / ray.d.y));
        deltaDist.z = length(mulS(ray.d, 1 / ray.d.z));
        const next = mul(add(maxV(vec3(0), step), add(mulS(mul(step, cf), -1), mul(step, c))), deltaDist);

        var ci = c.z * size * size + c.y * size + c.x;

        const closestHit = { index: -1, distance: Infinity };

        // Step through the grid while we're inside it
        while (!(c.x < 0 || c.y < 0 || c.z < 0 || c.x >= size || c.y >= size || c.z >= size)) {
            VoxelGrid.stepCount++;
            const vi = this._array[voxelsOff + ci];
            if (vi > 0) {
                const coff = childOff + this._array[childIndexOff + vi - 1];
                this.intersectLevel2(ray, coff, closestHit);
                if (closestHit.index >= 0) {
                    break;
                }
            }
            if (next.x < next.y) {
                if (next.x < next.z) {
                    next.x += deltaDist.x;
                    c.x += step.x;
                    ci += step.x;
                } else {
                    next.z += deltaDist.z;
                    c.z += step.z;
                    ci += step.z * size * size;
                }
            } else if (next.y < next.z) {
                next.y += deltaDist.y;
                c.y += step.y;
                ci += step.y * size;
            } else {
                next.z += deltaDist.z;
                c.z += step.z;
                ci += step.z * size * size;
            }
        }
        if (closestHit.index >= 0) {
            const hit = { obj: this._hit, distance: closestHit.distance };
            const index = closestHit.index;
            const off = 4 + index * 18;
            this._hit._e1.x = this._array[off + 3];
            this._hit._e1.y = this._array[off + 4];
            this._hit._e1.z = this._array[off + 5];
            this._hit._e2.x = this._array[off + 6];
            this._hit._e2.y = this._array[off + 7];
            this._hit._e2.z = this._array[off + 8];
            this._hit._normals[0].x = this._array[off + 9];
            this._hit._normals[0].y = this._array[off + 10];
            this._hit._normals[0].z = this._array[off + 11];
            this._hit._normals[1].x = this._array[off + 12];
            this._hit._normals[1].y = this._array[off + 13];
            this._hit._normals[1].z = this._array[off + 14];
            this._hit._normals[2].x = this._array[off + 15];
            this._hit._normals[2].y = this._array[off + 16];
            this._hit._normals[2].z = this._array[off + 17];
            return hit;
        }
        return null;
    }

    intersectLevel2(ray, coff, closestHit) {
        // Step through the voxel grid
        // On encountering a non-negative node,
        // look it up and intersect against it.

        // Do intersect against top-level VG bbox

        // Load _box from _array.
        // Add _box to _array.

        const headOff = coff;
        const origin = vec3(
            this._array[headOff],
            this._array[headOff + 1],
            this._array[headOff + 2]
        );
        const size = this._array[headOff + 3];
        const dims = vec3(
            this._array[headOff + 4],
            this._array[headOff + 5],
            this._array[headOff + 6]
        );
        const childSize = this._array[headOff + 7];
        const box = new Box(origin, add(origin, dims));
        const scale = dims.x / size;
        const voxelsOff = headOff + 8;
        const childOff = voxelsOff + size * size * size;

        const t = box.intersect(ray);
        if (t === null) {
            return null;
        }

        // Map hit to voxel coordinates
        const tray = new Ray(ray.o, ray.d, ray.time);
        tray.o = add(tray.o, mulS(tray.d, t.distance + scale * 0.0001));

        const cf = this.toGrid(tray.o, scale, origin);
        const c = floorV(cf);

        const deltaDist = vec3();
        const step = sign(ray.d);
        deltaDist.x = length(mulS(ray.d, 1 / ray.d.x));
        deltaDist.y = length(mulS(ray.d, 1 / ray.d.y));
        deltaDist.z = length(mulS(ray.d, 1 / ray.d.z));
        const next = mul(add(maxV(vec3(0), step), add(mulS(mul(step, cf), -1), mul(step, c))), deltaDist);

        var ci = c.z * size * size + c.y * size + c.x;

        // Step through the grid while we're inside it
        while (!(c.x < 0 || c.y < 0 || c.z < 0 || c.x >= size || c.y >= size || c.z >= size)) {
            VoxelGrid.stepCount++;
            const vi = this._array[voxelsOff + ci];
            if (vi > 0) {
                const coff = childOff + (vi - 1) * childSize;
                this.intersectTris(ray, coff, childSize, closestHit);
                if (closestHit.index >= 0) {
                    const px = floor((ray.o.x + ray.d.x * closestHit.distance - origin.x) / scale);
                    const py = floor((ray.o.y + ray.d.y * closestHit.distance - origin.y) / scale);
                    const pz = floor((ray.o.z + ray.d.z * closestHit.distance - origin.z) / scale);
                    if (px === c.x && py === c.y && pz === c.z) {
                        return;
                    }
                }
            }
            if (next.x < next.y) {
                if (next.x < next.z) {
                    next.x += deltaDist.x;
                    c.x += step.x;
                    ci += step.x;
                } else {
                    next.z += deltaDist.z;
                    c.z += step.z;
                    ci += step.z * size * size;
                }
            } else if (next.y < next.z) {
                next.y += deltaDist.y;
                c.y += step.y;
                ci += step.y * size;
            } else {
                next.z += deltaDist.z;
                c.z += step.z;
                ci += step.z * size * size;
            }
        }
    }


}

*/

void main() {
    gl_FragColor = vec4(0,0,0,0);    
}