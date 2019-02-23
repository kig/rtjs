#version 300 es

precision highp float;
precision highp int;

// WebGL 2 shader to ray trace.

uniform sampler2D arrayTex;
uniform float arrayTexWidth;

struct Hit {
    float index;
    float distance;
};

struct Ray {
    vec3 o;
    vec3 d;
    vec3 transmit;
    vec3 light;
    float pathLength;
    float lastTested;
};

struct Array {
    float width;
};

vec3 fromGrid(vec3 coord, float scale, vec3 origin) {
    return origin + (coord * scale);
}

vec3 toGrid(vec3 point, float scale, vec3 origin) {
    return (point - origin) / scale;
}

float readFloat(Array array, float index) {
    float v = floor(index / array.width);
    float u = index - (v * array.width);
    return texelFetch(arrayTex, ivec2(u, v), 0).r;
}

vec3 readVec3(Array array, float index) {
    return vec3(
        readFloat(array, index),
        readFloat(array, index + 1.0),
        readFloat(array, index + 2.0)
    );
}

float intersectBox(in Ray ray, in vec3 origin, in vec3 dims) {
    vec3 bmin = origin;
    vec3 bmax = origin + dims;
    vec3 invD = 1.0 / ray.d;
    vec3 t0 = (bmin - ray.o) * invD;
    vec3 t1 = (bmax - ray.o) * invD;
    vec3 swaps = (sign(invD) + 1.0) * 0.5;
    vec3 swaps2 = 1.0 - swaps;
    vec3 t0_ = t0 * swaps + t1 * swaps2;
    vec3 t1_ = t0 * swaps2 + t1 * swaps;

    float tmin = max(t0_.x, max(t0_.y, t0_.z));
    float tmax = min(t1_.x, max(t1_.y, t1_.z));

    if (tmax <= tmin) {
        return -1.0;
    }
    return tmin;
}

void intersectTri(in Array array, in Ray ray, in float triIndex, inout Hit closestHit) {
        float off = 4.0 + triIndex * 18.0;
        vec3 v0 = readVec3(array, off);
        vec3 e1 = readVec3(array, off+3.0);
        vec3 e2 = readVec3(array, off+6.0);

        vec3 h = cross(ray.d, e2);
        float a = dot(e1, h);

        if (a > -0.00000001 && a < 0.00000001) {
            return;
        }

        float f = 1.0 / a;
        vec3 s = ray.o - v0;
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


void intersectTris(in Array array, Ray ray, in float coff, in float childSize, inout Hit closestHit) {
    for (float j = 0.0; j < childSize; j++) {
        float triIndex = readFloat(array, coff + j);
        if (triIndex < 0.0) {
            break;
        }
        if (ray.lastTested != triIndex) {
            intersectTri(array, ray, triIndex, closestHit);
        }
    }
}

vec3 triNormal(in Array array, in vec3 point, in float triIndex) {
    float off = 4.0 + triIndex * 18.0;

    vec3 e1 = readVec3(array, off+3.0);
    vec3 e2 = readVec3(array, off+6.0);
    
    float u = dot(point, e1);
    float v = dot(point, e2);
    vec3 n0 = readVec3(array, off+9.0);
    vec3 n1 = readVec3(array, off+12.0);
    vec3 n2 = readVec3(array, off+15.0);

    return normalize(mix(mix(n0, n2, v), n1, u));
}


void intersectGridLeaf(in Array array, in Ray ray, in float headOff, inout Hit closestHit) {
    vec3 origin = readVec3(array, headOff);
    float size = readFloat(array, headOff + 3.0);
    vec3 dims = readVec3(array, headOff + 4.0);
    float childSize = readFloat(array, headOff + 7.0);
    float scale = dims.x / size;
    float voxelsOff = headOff + 8.0;
    float childOff = voxelsOff + size * size * size;

    float t = intersectBox(ray, origin, dims);
    if (t < 0.0) {
        return;
    }

    // Map hit to voxel coordinates
    Ray tray = ray;
    tray.o = tray.o + (tray.d * t + scale * 0.0001);

    vec3 cf = toGrid(tray.o, scale, origin);
    vec3 c = floor(cf);

    vec3 deltaDist = vec3(
        length(ray.d / ray.d.x),
        length(ray.d / ray.d.y),
        length(ray.d / ray.d.z)
    );
    vec3 cstep = sign(ray.d);
    vec3 next = (max(vec3(0.0), cstep) + cstep * (c - cf)) * deltaDist;

    float ci = c.z * size * size + c.y * size + c.x;

    // Step through the grid while we're inside it
    for (float i = 0.0; i < 3.0*size; i++) {
        if (c.x < 0.0 || c.y < 0.0 || c.z < 0.0 || c.x >= size || c.y >= size || c.z >= size) {
            return;
        }
        float vi = readFloat(array, voxelsOff + ci);
        if (vi > 0.0) {
            float coff = childOff + (vi - 1.0) * childSize;
            intersectTris(array, ray, coff, childSize, closestHit);
            if (closestHit.index >= 0.0) {
                vec3 p = toGrid(ray.o + ray.d * closestHit.distance, scale, origin);
                if (all(equal(p, c))) {
                    return;
                }
            }
        }
        if (next.x < next.y) {
            if (next.x < next.z) {
                next.x += deltaDist.x;
                c.x += cstep.x;
                ci += cstep.x;
            } else {
                next.z += deltaDist.z;
                c.z += cstep.z;
                ci += cstep.z * size * size;
            }
        } else if (next.y < next.z) {
            next.y += deltaDist.y;
            c.y += cstep.y;
            ci += cstep.y * size;
        } else {
            next.z += deltaDist.z;
            c.z += cstep.z;
            ci += cstep.z * size * size;
        }
    }
}


void intersectGridNode(in Array array, in Ray ray, inout Hit closestHit) {
    float headOff = 18.0 * readFloat(array, 0.0) + 4.0;

    vec3 origin = readVec3(array, headOff);
    float size = readFloat(array, headOff + 3.0);
    vec3 dims = readVec3(array, headOff + 4.0);
    float childSize = readFloat(array, headOff + 7.0);
    float scale = dims.x / size;
    float voxelsOff = headOff + 8.0;
    float childOff = voxelsOff + size * size * size;

    float t = intersectBox(ray, origin, dims);
    if (t < 0.0) {
        return;
    }

    // Map hit to voxel coordinates
    Ray tray = ray;
    tray.o = tray.o + (tray.d * t + scale * 0.0001);

    vec3 cf = toGrid(tray.o, scale, origin);
    vec3 c = floor(cf);

    vec3 deltaDist = vec3(
        length(ray.d / ray.d.x),
        length(ray.d / ray.d.y),
        length(ray.d / ray.d.z)
    );
    vec3 cstep = sign(ray.d);
    vec3 next = (max(vec3(0.0), cstep) + cstep * (c - cf)) * deltaDist;

    float ci = c.z * size * size + c.y * size + c.x;

    // Step through the grid while we're inside it
    for (float i = 0.0; i < 3.0*size; i++) {
        if (c.x < 0.0 || c.y < 0.0 || c.z < 0.0 || c.x >= size || c.y >= size || c.z >= size) {
            return;
        }
        float vi = readFloat(array, voxelsOff + ci);
        if (vi > 0.0) {
            float coff = childOff + (vi - 1.0) * childSize;
            intersectGridLeaf(array, ray, coff, closestHit);
            if (closestHit.index >= 0.0) {
                return;
            }
        }
        if (next.x < next.y) {
            if (next.x < next.z) {
                next.x += deltaDist.x;
                c.x += cstep.x;
                ci += cstep.x;
            } else {
                next.z += deltaDist.z;
                c.z += cstep.z;
                ci += cstep.z * size * size;
            }
        } else if (next.y < next.z) {
            next.y += deltaDist.y;
            c.y += cstep.y;
            ci += cstep.y * size;
        } else {
            next.z += deltaDist.z;
            c.z += cstep.z;
            ci += cstep.z * size * size;
        }
    }
}

