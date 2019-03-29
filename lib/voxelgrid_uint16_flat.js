
class Voxel3 {

    constructor() {
        this._objects = [];
        this._triangles = null;
        Voxel3.voxelCount++;
    }

    add(obj) {
        this._objects.push(obj);
        Voxel3.primitiveCount++;
        Voxel3.maxSize = max(Voxel3.maxSize, this._objects.length);
    }

    intersect(ray) {
        VoxelGrid.cmpCount += this._objects.length;
        let minT = Infinity;
        let hit = null;
        for (let i = 0; i < this._objects.length; i++) {
            const t = this._triangles[this._objects[i]].intersect(ray);
            if (t && t.distance < minT) {
                minT = t.distance;
                hit = t;
            }
        }
        return hit;
    }

}

Voxel3.primitiveCount = 0;
Voxel3.maxSize = 0;
Voxel3.voxelCount = 0;

class VoxelGrid3 {

    constructor(origin, dims, levelSizes, levelSizeIndex) {
        this._levelSize = levelSizes[levelSizeIndex];
        const size = this._levelSize;
        this._origin = origin;
        this._dims = dims;
        this._size = size;
        this._scale = this._dims.x / this._size;
        this._leaf = levelSizes.length - levelSizeIndex - 1;
        this._levelSizes = levelSizes;
        this._levelSizeIndex = levelSizeIndex;
        this._box = new Box(origin, add(origin, dims));
        this._voxels = new Int32Array(size * size * size);
        this._children = [];
    }

    serialize(root) {
        if (root === undefined) {
            root = true;
        }
        // const t0 = performance.now();
        let bigArray;
        const headerLength = 8;
        const bufs = [];
        if (this._leaf === 0) {
            const voxelSize = this._children.reduce((s, v) => max(v._objects.length, s), 0);
            const buf = new ArrayBuffer(headerLength * 4 + this._voxels.byteLength + 4 * voxelSize * this._children.length + 4);
            const header = new Float32Array(buf, 0);
            const iheader = new Int32Array(buf, 0);
            header[0] = this._origin.x;
            header[1] = this._origin.y;
            header[2] = this._origin.z;
            iheader[3] = this._levelSize;
            header[4] = this._dims.x;
            header[5] = this._dims.y;
            header[6] = this._dims.z;
            iheader[7] = voxelSize;
            let i32 = new Int32Array(buf, headerLength * 4);
            for (let i = 0; i < i32.length; i++) {
                i32[i] = -1;
            }
            for (let i = 0; i < this._voxels.length; i++) {
                const v = this._voxels[i];
                if (v !== 0) {
                    if (this._children[v-1]._objects.length === 1) {
                        this._voxels[i] = -(this._children[v-1]._objects[0] + 1);
                    }
                }
            }
            i32.set(this._voxels);
            const triArray = new Int32Array(buf, headerLength * 4 + this._voxels.byteLength, this._children.length * voxelSize + 1);
            // triArray[0] = this._children.length;
            // console.log(this._children.length, voxelSize);
            // const lengths = [];
            for (let i = 0; i < this._children.length; i++) {
                // lengths.push(this._children[i]._objects.length);
                triArray.set(this._children[i]._objects, voxelSize * i);
            }
            // console.log(lengths.join(", "));
            bufs.push(header);
        } else {
            const buf = new ArrayBuffer(headerLength * 4 + this._voxels.byteLength + 4 * this._children.length);
            const header = new Float32Array(buf, 0);
            const iheader = new Int32Array(buf, 0);
            header[0] = this._origin.x;
            header[1] = this._origin.y;
            header[2] = this._origin.z;
            iheader[3] = -this._levelSize;
            header[4] = this._dims.x;
            header[5] = this._dims.y;
            header[6] = this._dims.z;
            iheader[7] = this._children.length;
            let i32 = new Int32Array(buf, headerLength * 4);
            i32.set(this._voxels);
            const cArray = new Int32Array(buf, headerLength * 4 + this._voxels.byteLength, this._children.length);
            bufs.push(header);
            let childOffset = 0;
            for (let i = 0; i < this._children.length; i++) {
                cArray[i] = childOffset;
                // console.log(childOffset);
                const cbuf = this._children[i].serialize(false);
                childOffset += cbuf.length;
                bufs.push(cbuf);
            }
        }
        if (root) {
            const triArr = new Float32Array(4 + this._triangles.length * 3 * (3 + 3)); // 3 verts, 3 normals
            const itriArr =  new Int32Array(triArr.buffer);
            itriArr[0] = this._triangles.length;
            for (let i = 0, j = 4; i < this._triangles.length; ++i, ++j) {
                const tri = this._triangles[i];
                const v = tri._vertices[0];
                const e1 = tri._e1;
                const e2 = tri._e2;
                const n0 = tri._normals[0];
                const n1 = tri._normals[1];
                const n2 = tri._normals[2];
                triArr[j] = v.x;
                triArr[++j] = v.y;
                triArr[++j] = v.z;
                triArr[++j] = e1.x;
                triArr[++j] = e1.y;
                triArr[++j] = e1.z;
                triArr[++j] = e2.x;
                triArr[++j] = e2.y;
                triArr[++j] = e2.z;
                triArr[++j] = n0.x;
                triArr[++j] = n0.y;
                triArr[++j] = n0.z;
                triArr[++j] = n1.x;
                triArr[++j] = n1.y;
                triArr[++j] = n1.z;
                triArr[++j] = n2.x;
                triArr[++j] = n2.y;
                triArr[++j] = n2.z;
            }
            bufs.unshift(triArr);
        }
        // console.log(performance.now() - t0);
        bigArray = new Float32Array(bufs.reduce((s, i) => s + i.length, 0));
        bufs.reduce((offset, arr) => {
            bigArray.set(arr, offset);
            return offset + arr.length;
        }, 0);
        return bigArray;
    }

    add(index) {
        const obj = this._triangles[index];
        const size = this._size;
        const scaleVec = vec3(this._scale);
        const voxelVec = vec3(0,0,0);
        const bbox = obj.getBoundingBox();
        bbox.mid = add(bbox.min, mulS(scaleVec, 0.5));
        const start = maxV(vec3(0), minV(vec3(this._size - 1), floorV(this.toGrid(bbox.min))));
        const end = maxV(vec3(0), minV(vec3(this._size - 1), floorV(this.toGrid(bbox.max))));
        for (let z = start.z; z <= end.z; z++) {
            const cz = z * size * size;
            for (let y = start.y; y <= end.y; y++) {
                const cy = cz + y * size;
                for (let x = start.x; x <= end.x; x++) {
                    const c = cy + x;
                    voxelVec.x = x;
                    voxelVec.y = y;
                    voxelVec.z = z;
                    bbox.min = this.fromGrid(voxelVec);
                    bbox.max.x = bbox.min.x + this._scale;
                    bbox.max.y = bbox.min.y + this._scale;
                    bbox.max.z = bbox.min.z + this._scale;
                    bbox.mid.x = bbox.min.x + this._scale * 0.5;
                    bbox.mid.y = bbox.min.y + this._scale * 0.5;
                    bbox.mid.z = bbox.min.z + this._scale * 0.5;
                    let intersects = false;
                    // if (this._leaf === 0) {
                        intersects = obj.intersectBoxFast(bbox) && obj.intersectBox(bbox);
                    // } else {
                        // intersects = obj.intersectBoxFast(bbox);
                    // }
                    if (intersects) {
                        let vi = this._voxels[c];
                        if (!vi) {
                            vi = this._voxels[c] = this._children.length + 1;
                            var v = this._leaf === 0 ? new Voxel3() : new VoxelGrid3(bbox.min, scaleVec, this._levelSizes, this._levelSizeIndex + 1);
                            this._children.push(v);
                        }
                        var v = this._children[vi - 1];
                        v._triangles = this._triangles;
                        v.add(index);
                    }
                }
            }
        }
    }

    addTriangles(tris) {
        this._triangles = tris.slice();
        for (var i = 0; i < tris.length; i++) {
            this.add(i);
        }
    }
    
    prune() {
        this._children = this._children.map(c => c.prune ? c.prune() : c);
        if (this._children.every(c => c._objects && c._objects.length === 1 && c._objects[0] === this._children[0]._objects[0])) {
            VoxelGrid2.pruned += this._children.length;
            return this._children[0];	
        }
        return this;
    }
    
    fromGrid(coord) {
        return add(this._origin, mulS(coord, this._scale));
    }

    toGrid(point) {
        return divS(sub(point, this._origin), this._scale);
    }

    intersect(ray) {
        const t = this._box.intersect(ray);
        if (t === null) {
            return null;
        }

        const tray = new Ray(ray.o, ray.d, ray.time);
        tray.o = add(tray.o, mulS(tray.d, t.distance + this._scale * 0.0001));

        const cf = this.toGrid(tray.o);
        const c = floorV(cf);

        const deltaDist = vec3();
        const step = sign(ray.d);
        deltaDist.x = length(mulS(ray.d, 1 / ray.d.x));
        deltaDist.y = length(mulS(ray.d, 1 / ray.d.y));
        deltaDist.z = length(mulS(ray.d, 1 / ray.d.z));
        const next = mul(add(maxV(vec3(0), step), add(mulS(mul(step, cf), -1), mul(step, c))), deltaDist);
        const size = this._size;
        var ci = c.z * size * size + c.y * size + c.x;
        while (!(c.x < 0 || c.y < 0 || c.z < 0 || c.x >= size || c.y >= size || c.z >= size)) {
            VoxelGrid.stepCount++;
            const vi = this._voxels[ci];
            if (vi) {
                const voxel = this._children[vi - 1];
                const hit = voxel.intersect(ray);
                if (hit) {
                    const px = floor((ray.o.x + ray.d.x * hit.distance - this._origin.x) / this._scale);
                    const py = floor((ray.o.y + ray.d.y * hit.distance - this._origin.y) / this._scale);
                    const pz = floor((ray.o.z + ray.d.z * hit.distance - this._origin.z) / this._scale);
                    if (px === c.x && py === c.y && pz === c.z) {
                        return hit;
                    }
                    // const p = add(ray.o, mulS(ray.d, hit.distance));
                    // if (all(equal(c, floorV(this.toGrid(p))))) {
                    // 	return hit;
                    // }
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
        return null;
    }
}



/*
    Tracing from a flat Float32Array.
    This is to prepare for a WebGL port.

    With WebGL the main challenge might be how to do the grid stepping with minimal divergence.
    Especially since the voxel grid has two layers, and stepping down to the lower layer
    screws up rays that don't.
    
*/
class SerializedVG {
    constructor(vgFloat32Array, color, isLeaf) {
        this._color = color;
        this._array = vgFloat32Array;
        this._iarray = new Int32Array(vgFloat32Array.buffer);
        this._isLeaf = isLeaf;
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

    fromGrid(coord, scale, origin) {
        return add(origin, mulS(coord, scale));
    }

    toGrid(point, scale, origin) {
        return divS(sub(point, origin), scale);
    }

    intersectTri(ray, triIndex, closestHit) {
        const off = 4 + triIndex * 18;
        const v0 = vec3(
            this._array[off],
            this._array[off + 1],
            this._array[off + 2]
        );
        const e1 = vec3(
            this._array[off + 3],
            this._array[off + 4],
            this._array[off + 5]
        );
        const e2 = vec3(
            this._array[off + 6],
            this._array[off + 7],
            this._array[off + 8]
        );

        const h = cross(ray.d, e2);
        const a = dot(e1, h);

        if (a > -0.00000001 && a < 0.00000001) {
            return;
        }

        const f = 1.0 / a;
        const s = sub(ray.o, v0);
        const u = f * dot(s, h);

        if (u < -0.00000001 || u > 1.00000001) {
            return;
        }

        const q = cross(s, e1);
        const v = f * dot(ray.d, q);

        if (v < -0.00000001 || u + v > 1.00000002) {
            return;
        }

        // at this stage we can compute t to find out where
        // the intersection point is on the line
        const t = f * dot(e2, q);

        if (t <= 0.00000001 || t >= closestHit.distance) {
            return;
        }

        closestHit.index = triIndex;
        closestHit.distance = t;
    }

    intersectTris(ray, coff, childSize, closestHit) {
        for (let j = coff; j < coff + childSize; j++) {
            const triIndex = this._iarray[j];
            if (triIndex < 0) {
                break;
            }
            VoxelGrid.cmpCount++;
            this.intersectTri(ray, triIndex, closestHit);
        }
        return closestHit;
    }

    intersect(ray) {
        const closestHit = { index: -1, distance: Infinity };

        this.intersectGridNode(ray, this._iarray[0] * 18 + 4, closestHit)

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

    intersectGridNode(ray, headOff, closestHit) {
        const origin = vec3(
            this._array[headOff],
            this._array[headOff + 1],
            this._array[headOff + 2]
        );
        const osize = this._iarray[headOff + 3];
        const isLeaf = osize > 0;
        const size = abs(osize);
        const dims = vec3(
            this._array[headOff + 4],
            this._array[headOff + 5],
            this._array[headOff + 6]
        );
        const childCount = this._iarray[headOff + 7];
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

        const deltaDist = vec3(
            length(mulS(ray.d, 1 / ray.d.x)),
            length(mulS(ray.d, 1 / ray.d.y)),
            length(mulS(ray.d, 1 / ray.d.z))
        );
        const step = sign(ray.d);
        const next = mul( deltaDist, add(maxV(vec3(0), step), sub(mul(step, c), mul(step, cf))) );

        var ci = c.z * size * size + c.y * size + c.x;

        // Step through the grid while we're inside it
        while (!(c.x < 0 || c.y < 0 || c.z < 0 || c.x >= size || c.y >= size || c.z >= size)) {
            VoxelGrid.stepCount++;
            const vi = this._iarray[voxelsOff + ci];
            if (vi > 0) {
                if (isLeaf) {
                    const coff = childIndexOff + (vi - 1) * childCount;
                    this.intersectTris(ray, coff, childCount, closestHit);
                    if (closestHit.index >= 0) {
                        const px = floor((ray.o.x + ray.d.x * closestHit.distance - origin.x) / scale);
                        const py = floor((ray.o.y + ray.d.y * closestHit.distance - origin.y) / scale);
                        const pz = floor((ray.o.z + ray.d.z * closestHit.distance - origin.z) / scale);
                        if (px === c.x && py === c.y && pz === c.z) {
                            return;
                        }
                    }
                } else {
                    const coff = childOff + this._iarray[childIndexOff + vi - 1];
                    this.intersectGridNode(ray, coff, closestHit);
                    if (closestHit.index >= 0) {
                        break;
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
