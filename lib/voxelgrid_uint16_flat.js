
class Voxel3 {

    constructor() {
        this._objects = [];
        this._triangles = null;
    }

    add(obj) {
        this._objects.push(obj);
        Voxel.primitiveCount++;
        Voxel.maxSize = max(Voxel.maxSize, this._objects.length);
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
        this._voxels = new Uint32Array(size * size * size);
        this._children = [];
    }

    serialize(root) {
        if (root === undefined) {
            root = true;
        }
        // const t0 = performance.now();
        let blob;
        if (this._leaf === 0) {
            const voxelSize = this._children.reduce((s, v) => max(v._objects.length, s), 0);
            const buf = new ArrayBuffer(12 + this._voxels.byteLength + 4 * voxelSize * this._children.length);
            const header = new Float32Array(buf, 0, 3);
            header[0] = this._levelSize;
            header[1] = this._children.length;
            header[2] = voxelSize;
            new Float32Array(buf, 8).set(this._voxels);
            const triArray = new Float32Array(buf, this._voxels.byteLength + 12, this._children.length * voxelSize);
            // console.log(this._children.length, voxelSize);
            // const lengths = [];
            for (let i = 0; i < this._children.length; i++) {
                // lengths.push(this._children[i]._objects.length);
                triArray.set(this._children[i]._objects, voxelSize * i);
            }
            // console.log(lengths.join(", "));
            blob = new Blob([buf]);
        } else {
            const buf = new ArrayBuffer(12 + this._voxels.byteLength + 4 * this._children.length);
            const header = new Float32Array(buf, 0, 3);
            header[0] = this._levelSize;
            header[1] = this._children.length;
            header[2] = this._children[0]._size;
            new Float32Array(buf, 12).set(this._voxels);
            const cArray = new Float32Array(buf, this._voxels.byteLength + 12, this._children.length);
            const bufs = [buf];
            if (root) {
                const triArr = new Float32Array(4 + this._triangles.length * 3 * (3+3)); // 3 verts, 3 normals
                triArr[0] = this._triangles.length;
                for (let i = 0, j = 0; i < this._triangles.length; i++, j++) {
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
                bufs.push(triArr.buffer);
            }
            let childOffset = 0;
            for (let i = 0; i < this._children.length; i++) {
                cArray[i] = childOffset;
                const cbuf = this._children[i].serialize(false);
                childOffset += cbuf.size;
                bufs.push(cbuf);
            }
            blob = new Blob(bufs);
        }
        if (this._root) {

        }
        // console.log(performance.now() - t0);
        return blob;
    }

    add(index) {
        const obj = this._triangles[index];
        const size = this._size;
        const bbox = obj.getBoundingBox();
        const start = maxV(vec3(0), minV(vec3(this._size - 1), floorV(this.toGrid(bbox.min))));
        const end = maxV(vec3(0), minV(vec3(this._size - 1), floorV(this.toGrid(bbox.max))));
        for (let z = start.z; z <= end.z; z++) {
            const cz = z * size * size;
            for (let y = start.y; y <= end.y; y++) {
                const cy = cz + y * size;
                for (let x = start.x; x <= end.x; x++) {
                    const c = cy + x;
                    let vi = this._voxels[c];
                    if (!vi) {
                        vi = this._voxels[c] = this._children.length + 1;
                        var v = this._leaf === 0 ? new Voxel3() : new VoxelGrid3(this.fromGrid(vec3(x, y, z)), vec3(this._scale), this._levelSizes, this._levelSizeIndex + 1);
                        this._children.push(v);
                    }
                    var v = this._children[vi - 1];
                    v._triangles = this._triangles;
                    v.add(index);
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
