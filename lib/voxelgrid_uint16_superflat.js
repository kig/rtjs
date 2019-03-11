var VoxelGrid4 = function (origin, dims, levelSizes, levelSizeIndex, voxelIndex, viOffset, voxels, voOffset) {
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
	this._voxelIndex = voxelIndex || new Uint32Array(
        (((1 + size * size * size) * size * size * size) >> 5) +
        (size * size * size) * size * size * size
    );
	this._voxels = voxels || new Uint32Array(this._voxelIndex.buffer, (((1 + size * size * size) * size * size * size) >> 5)*4);
	this._voxelIndexOffset = viOffset || 0;
	this._voxelsOffset = voOffset || 0;
	this._children = [];
	VoxelGrid.count++;
};

VoxelGrid.count = 0;
VoxelGrid4.cacheLineLoadCount = 0;
VoxelGrid4.cacheLines = {};

VoxelGrid4.prototype.add = function (obj) {
	var size = this._size;
	var bbox = obj.getBoundingBox();
	var scaleVec = vec3(this._scale);
	var start = maxV(vec3(0), minV(vec3(this._size - 1), floorV(this.toGrid(this._origin, bbox.min, this._scale)))); //sub(bbox.min, mulS(this._dims, (this._levelSizeIndex===0?1:0) * 1.01/size))))));
	var end = maxV(vec3(0), minV(vec3(this._size - 1), floorV(this.toGrid(this._origin, bbox.max, this._scale)))); //(add(bbox.max, mulS(this._dims, (this._levelSizeIndex===0?1:0) * 1.01/size))))));
	for (var z = start.z; z <= end.z; z++) {
		var cz = z * size * size;
		for (var y = start.y; y <= end.y; y++) {
			var cy = cz + y * size;
			for (var x = start.x; x <= end.x; x++) {
				var c = cy + x;
				bbox.min = this.fromGrid(this._origin, vec3(x,y,z), this._scale);
				bbox.max = add(bbox.min, scaleVec);
				if (obj.intersectBox(bbox)) {
					var vi = (this._voxelIndex[this._voxelIndexOffset + (c >> 5)] >> (c & 31)) & 1;
					if (!vi) {
						this._voxelIndex[this._voxelIndexOffset + (c >> 5)] |= (1 << (c & 31));
						if (this._leaf === 0) {
							this._voxels[this._voxelsOffset + c] = Voxel.all.length;
							new Voxel();
						}
					}
					if (this._leaf === 0) {
                        Voxel.all[this._voxels[this._voxelsOffset + c]].add(obj);
					} else {
						(new VoxelGrid4(
                            this.fromGrid(this._origin, vec3(x, y, z), this._scale), 
                            scaleVec, 
                            this._levelSizes, 
                            this._levelSizeIndex + 1,
                            this._voxelIndex,
                            this._voxelIndexOffset + (((1+c) * size*size*size) >> 5),
                            this._voxels,
                            this._voxelsOffset + (c * size*size*size)
                        )).add(obj);
					}
				}
			}
		}
	}
};

VoxelGrid4.prototype.prune = function() {
	this._children = this._children.map(c => c.prune ? c.prune() : c);
	if (this._children.every(c => c._objects && c._objects.length === 1 && c._objects[0] === this._children[0]._objects[0])) {
		VoxelGrid2.pruned += this._children.length;
		return this._children[0];	
	}
	return this;
};

VoxelGrid4.prototype.fromGrid = function (origin, coord, scale) {
	return add(origin, mulS(coord, scale));
};

VoxelGrid4.prototype.toGrid = function (origin, point, scale) {
	return divS(sub(point, origin), scale);
};

VoxelGrid4.prototype.getSubNode = function (c, ci, size, ray, step, deltaDist) {
    const voxelIndexOffset = this._voxelIndexOffset + (((1+ci) * size*size*size) >> 5);
    const voxelsOffset = this._voxelsOffset + (ci * size*size*size);

    const childOrigin = add(this._origin, mulS(c, this._scale));

    const childBox = new Box(
        childOrigin,
        addS(childOrigin, this._scale)
    );
	const t = childBox.intersect(ray);
	if (t === null) {
		return null;
	}

	const cf = this.toGrid(childOrigin, add(ray.o, mulS(ray.d, t.distance + this._scale / size * 0.0001)), this._scale / size);
	const cc = floorV(cf);

	const cci = cc.z * size * size + cc.y * size + cc.x;

    const next = mul(add(maxV(vec3(0), step), add(mulS(mul(step, cf), -1), mul(step, cc))), deltaDist);

    return [cc, next, cci, voxelIndexOffset, voxelsOffset];
};

VoxelGrid4.prototype.intersect = function (ray) {
    VoxelGrid.descendCount++;

	const t = this._box.intersect(ray);
	if (t === null) {
		return null;
	}

    const deltaDist = vec3();
	const step = sign(ray.d);
	deltaDist.x = length(mulS(ray.d, 1 / ray.d.x));
	deltaDist.y = length(mulS(ray.d, 1 / ray.d.y));
	deltaDist.z = length(mulS(ray.d, 1 / ray.d.z));

	const cf = this.toGrid(this._origin, add(ray.o, mulS(ray.d, t.distance + this._scale * 0.0001)), this._scale);
	var c = floorV(cf);

	const size = this._size;
	var ci = c.z * size * size + c.y * size + c.x;

    var next = mul(add(maxV(vec3(0), step), add(mulS(mul(step, cf), -1), mul(step, c))), deltaDist);
   
    var voxelIndexOffset = this._voxelIndexOffset;
    var voxelsOffset = this._voxelsOffset;

    var topLevel = [c, next, ci, voxelIndexOffset, voxelsOffset];
    var onTopLevel = true;

	while (true) {
        if (c.x < 0 || c.y < 0 || c.z < 0 || c.x >= size || c.y >= size || c.z >= size) {
            if (onTopLevel) {
                // escaped top level
                break;
            } else {
                // pop stack
                onTopLevel = true;
                [c, next, ci, voxelIndexOffset, voxelsOffset] = topLevel;
            }
        } else {
            VoxelGrid.stepCount++;
            ray.light.y += 0.02;
            if (!VoxelGrid4.cacheLines[(voxelIndexOffset + (ci >> 5)) >> 3]) {
                VoxelGrid4.cacheLines[(voxelIndexOffset + (ci >> 5)) >> 3] = 1;
                VoxelGrid4.cacheLineLoadCount++;
            }
            const vi = (this._voxelIndex[voxelIndexOffset + (ci >> 5)] >> (ci & 31)) & 1;
            if (vi !== 0) {
                // When at top-level node, jump to second-level node and continue intersect
                // (push current node on stack, replace with sub-node stuff)
                if (onTopLevel) {
                    VoxelGrid.descendCount++;
                    topLevel = [c, next, ci, voxelIndexOffset, voxelsOffset];
                    [c, next, ci, voxelIndexOffset, voxelsOffset] = this.getSubNode(c, ci, size, ray, step, deltaDist);
                    onTopLevel = false;
                    continue;
                } else {
                    // When at second-level node, jump to voxel
                    // if (!VoxelGrid4.cacheLines[((((1 + size * size * size) * size * size * size) >> 5) + voxelsOffset + ci) >> 3]) {
                    //     VoxelGrid4.cacheLines[((((1 + size * size * size) * size * size * size) >> 5) + voxelsOffset + ci) >> 3] = 1;
                    //     VoxelGrid4.cacheLineLoadCount++;
                    // }
                    if (!VoxelGrid4.cacheLines[((((1 + size * size * size) * size * size * size) >> 5) + size**6 + this._voxels[voxelsOffset + ci]*8) >> 3]) {
                        VoxelGrid4.cacheLines[((((1 + size * size * size) * size * size * size) >> 5) + size**6 + this._voxels[voxelsOffset + ci]*8) >> 3] = 1;
                        VoxelGrid4.cacheLineLoadCount++;
                    }
                    const voxel = Voxel.all[this._voxels[voxelsOffset + ci]];
                    const hit = voxel.intersect(ray);
                    // VoxelGrid4.cacheLineLoadCount += voxel._objects.length;
                    if (hit) {
                        // Can do this test by comparing hit distance to bbox intersect tmin & tmax too
                        return hit;
                        const px = floor((ray.o.x + ray.d.x * hit.distance - this._origin.x) / this._scale);
                        const py = floor((ray.o.y + ray.d.y * hit.distance - this._origin.y) / this._scale);
                        const pz = floor((ray.o.z + ray.d.z * hit.distance - this._origin.z) / this._scale);
                        if (px === c.x && py === c.y && pz === c.z) {
                            return hit;
                        }
                    }
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
	return null;
};
