var VoxelGrid4 = function (origin, dims, levelSizes, levelSizeIndex, voxelIndex, viOffset, voxels, voOffset, triIndices, triangles) {
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
	this._voxelIndex = voxelIndex || new Int32Array(
        (((1 + size * size * size) * size * size * size) >> 5)
    );
	this._voxels = voxels || new Uint32Array((size ** 6) * 16);
	this._voxelIndexOffset = viOffset || 0;
	this._voxelsOffset = voOffset || 0;
	this._triIndices = triIndices || [];

	// Instead of storing the entire grid, store just full nodes.
	// To calculate memory position of a node, count the bits set before it.
	// For each grid plane, store a prefix sum of set bits.
	//
	// This is best done as a post-processing stage:
	//  1) Generate a sparse grid with JS objects.
	//  2) Compress it into dense lists of grids and voxels.
	//  3) Keep a rolling counter of grid memory locations and voxel memory locations.
	//  4) Store counter after each plane.
	//
	// To traverse the grid:
	//  1) Do your normal traversal to find a non-zero grid coordinate (x,y,z).
	//  2) Fetch the prefix sum for the coordinate plane z.
	//  3) Count bits set on plane z before (x,y,z), add to prefix sum.
	//  4) Jump to subgrid / triangle list address calculated from the bit count.

	this._triangles = triangles;
	VoxelGrid4.count++;
};

VoxelGrid4.count = 0;
VoxelGrid4.cacheLineLoadCount = 0;
VoxelGrid4.cacheLines = {};
VoxelGrid4.cacheLineGroups = {};
VoxelGrid4.memoryAccessesGroups = {};
VoxelGrid4.memoryAccesses = 0;

VoxelGrid4.memoryAccess = function(addr, bytes, name) {
	const startCL = (addr) >> 5;
	const endCL = (addr + bytes-1) >> 5;
	for (let i = startCL; i <= endCL; i++) {
		if (!VoxelGrid4.cacheLines[i]) {
			VoxelGrid4.cacheLines[i] = 1;
			VoxelGrid4.cacheLineLoadCount++;
			VoxelGrid4.cacheLineGroups[name] = (VoxelGrid4.cacheLineGroups[name] || 0) + 1;
		}
	}
	VoxelGrid4.memoryAccesses += bytes;
	VoxelGrid4.memoryAccessesGroups[name] = (VoxelGrid4.memoryAccessesGroups[name] || 0) + bytes;
};

VoxelGrid4.prototype.serialize = function() { 
	const triArr = new Float32Array(this._triangles.length * 9);
	const nmlArr = new Float32Array(this._triangles.length * 6);
	for (var i = 0, j = 0, k = 0; i < this._triangles.length; i++, j++, k++) {
		const tri = this._triangles[i];
		const v0 = tri._vertices[0];
		const v1 = tri._vertices[1];
		const v2 = tri._vertices[2];
		const n0 = tri._normals[0];
		const n1 = tri._normals[1];
		const n2 = tri._normals[2];

		triArr[j] = v0.x;
		triArr[++j] = v0.y;
		triArr[++j] = v0.z;
		triArr[++j] = v1.x;
		triArr[++j] = v1.y;
		triArr[++j] = v1.z;
		triArr[++j] = v2.x;
		triArr[++j] = v2.y;
		triArr[++j] = v2.z;

		let uv = toOct(n0); 
		nmlArr[k] = uv.x;
		nmlArr[++k] = uv.y;
		uv = toOct(n1);
		nmlArr[++k] = uv.x;
		nmlArr[++k] = uv.y;
		uv = toOct(n2); 
		nmlArr[++k] = uv.x;
		nmlArr[++k] = uv.y;
	}
	this._triArr = triArr;
	this._nmlArr = nmlArr;
	return {
		voxelIndex: this._voxelIndex,
		triIndices: this._voxels,
		// triIndices: this._triIndicesFlat,
		triangles: triArr,
		normals: nmlArr
	};
};

VoxelGrid4.prototype.addTriangles = function(tris) {
	this._triangles = tris.slice();
	for (let i = 0; i < tris.length; i++) {
		this.add(i);
	}
	this._triIndicesFlat = new Uint16Array(16 * this._triIndices.length);
	for (let i = 0; i < this._triIndices.length; i++) {
		var idxs = this._triIndices[i];
		for (let j = 0; j < idxs.length; j++) {
			this._triIndicesFlat[i*16 + j] = idxs[j] + 1;
		}
	}
};

VoxelGrid4.prototype.add = function (index) {
	const obj = this._triangles[index];
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
				if (obj.intersectBoxFast(bbox) && obj.intersectBox(bbox)) {
					var vi = (this._voxelIndex[this._voxelIndexOffset + (c >> 5)] >> (c & 31)) & 1;
					if (!vi) {
						this._voxelIndex[this._voxelIndexOffset + (c >> 5)] |= (1 << (c & 31));
						// if (this._leaf === 0) {
						// 	this._voxels[this._voxelsOffset + c] = 0;
						// }
					}
					if (this._leaf === 0) {
						let off = this._voxelsOffset + c*16;
						while (this._voxels[off] > 0) {
							++off;
						}
						this._voxels[off] = index+1;
						// if (vx === 0) {
						// 	this._voxels[this._voxelsOffset + c] = index | 0x80000000;
						// } else if (vx & 0x80000000 && !(vx & 0x40000000)) {
						// 	this._voxels[this._voxelsOffset + c] = vx | (index << 15) | 0x40000000;
						// } else {
						// 	if (vx & 0x80000000) {
						// 		this._voxels[this._voxelsOffset + c] = this._triIndices.length;
						// 		if (vx & 0x40000000) {
						// 			this._triIndices.push([vx & 0x7fff, (vx >> 15) & 0x7fff]);
						// 		} else {
						// 			this._triIndices.push([vx & ~0x80000000]);
						// 		}
						// 		vx = this._triIndices.length-1;
						// 	}
						// 	this._triIndices[vx].push(index);
						// }
					} else {
						(new VoxelGrid4(
                            this.fromGrid(this._origin, vec3(x, y, z), this._scale), 
                            scaleVec, 
                            this._levelSizes, 
                            this._levelSizeIndex + 1,
                            this._voxelIndex,
                            this._voxelIndexOffset + (((1+c) * size*size*size) >> 5),
                            this._voxels,
							this._voxelsOffset + (c * size*size*size * 16),
							this._triIndices,
							this._triangles
                        )).add(index);
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
    const voxelsOffset = this._voxelsOffset + (ci * size*size*size) * 16;

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

    return [cc, next, cci, voxelIndexOffset, voxelsOffset, childOrigin, this._scale / size];
};

VoxelGrid4.prototype.triIntersect = function(index, ray) {
	VoxelGrid4.memoryAccess(1e7 + (index*9*4), 9*4, 'tri');
	return this._triangles[index].intersect(ray);
	const off = index * 9;
	const v0 = vec3(
		this._triArr[off],
		this._triArr[off + 1],
		this._triArr[off + 2]
	);
	const e1 = sub(vec3(
		this._triArr[off + 3],
		this._triArr[off + 4],
		this._triArr[off + 5]
	), v0);
	const e2 = sub(vec3(
		this._triArr[off + 6],
		this._triArr[off + 7],
		this._triArr[off + 8]
	), v0);

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

	if (t <= 0.00000001) {
		return;
	}

	return {index: index, distance: t, obj: this._triangles[index]};
};

VoxelGrid4.prototype.triListIntersect = function (indexArray, indexCoord, ray) {
	var minT = Infinity;
	var hit = null;
	// Voxel.visitCount++;
	let i = 0;
	for (i = 0; i < 16; i++) {
		VoxelGrid4.memoryAccess(2e7 + ((indexCoord+i)*2), 2, 'triList');
		let index = indexArray[indexCoord+i];
		if (index === 0) {
			break;
		}
		index--;
		ray.light.x += 0.1;
		// VoxelGrid.cmpCount++;
		const t = this.triIntersect(index, ray);
		if (t && t.distance < minT) {
			minT = t.distance;
			hit = t;
		}
	}
	// Voxel.visitCounts[i] = (Voxel.visitCounts[i] || 0) + 1;
	return hit;
};

VoxelGrid4.prototype.intersect = function (ray) {
    // VoxelGrid.descendCount++;

	VoxelGrid4.memoryAccess(3e7, 8*4, 'header');

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
	var origin = this._origin;
	var scale = this._scale;

    var topLevel = [c, next, ci, voxelIndexOffset, voxelsOffset, origin, scale];
    var onTopLevel = true;

	while (true) {
        if (c.x < 0 || c.y < 0 || c.z < 0 || c.x >= size || c.y >= size || c.z >= size) {
            if (onTopLevel) {
                // escaped top level
                break;
            } else {
                // pop stack
                onTopLevel = true;
                [c, next, ci, voxelIndexOffset, voxelsOffset, origin, scale] = topLevel;
            }
        } else {
            // VoxelGrid.stepCount++;
            ray.light.y += 0.02;
			VoxelGrid4.memoryAccess(4*(voxelIndexOffset + (ci >> 5)), 4, 'vg');
            const vi = (this._voxelIndex[voxelIndexOffset + (ci >> 5)] >> (ci & 31)) & 1;
            if (vi !== 0) {
                // When at top-level node, jump to second-level node and continue intersect
                // (push current node on stack, replace with sub-node stuff)
                if (onTopLevel) {
                    // VoxelGrid.descendCount++;
                    topLevel = [c, next, ci, voxelIndexOffset, voxelsOffset, origin, scale];
                    [c, next, ci, voxelIndexOffset, voxelsOffset, origin, scale] = this.getSubNode(c, ci, size, ray, step, deltaDist);
                    onTopLevel = false;
                    continue;
                } else {
                    // When at second-level node, jump to voxel
					// VoxelGrid4.memoryAccess(4*((((1 + size**3) * size**3) >> 5) + voxelsOffset + ci), 4, 'voxel');
					// const voxel = this._voxels[voxelsOffset + ci];
					// let hit;
					// if (voxel & 0x80000000) {
					// 	VoxelGrid.cmpCount++;
					// 	ray.light.y += 0.125;
					// 	ray.light.z += 0.25;
					// 	hit = this.triIntersect(voxel & 0x7fff, ray);
					// 	if (voxel & 0x40000000) {
					// 		VoxelGrid.cmpCount++;
					// 		ray.light.y += 0.25;
					// 		ray.light.z += 0.25;
					// 		let hit2 = this.triIntersect((voxel >> 15) & 0x7fff, ray);
					// 		if (hit2 && (!hit || hit.distance > hit2.distance)) {
					// 			hit = hit2;
					// 		}
					// 	}
					// } else {
					// 	hit = this.triListIntersect(this._triIndicesFlat, voxel * 16, ray);
					// }
					const hit = this.triListIntersect(this._voxels, voxelsOffset + ci * 16, ray);
					if (hit) {
                        // Can do this test by comparing hit distance to bbox intersect tmin & tmax too
                        const px = floor((ray.o.x + ray.d.x * hit.distance - origin.x) / scale);
                        const py = floor((ray.o.y + ray.d.y * hit.distance - origin.y) / scale);
                        const pz = floor((ray.o.z + ray.d.z * hit.distance - origin.z) / scale);
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
