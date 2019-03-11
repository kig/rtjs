var VoxelGrid2 = function (origin, dims, levelSizes, levelSizeIndex, voxelIndex, viOffset, voxels, voOffset) {
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
	this._voxelIndex = voxelIndex || new Uint32Array(Math.ceil((1 + size * size * size) * size * size * size / 32));
	this._voxels = voxels || new Uint32Array((1 + size * size * size) * size * size * size);
	this._voxelIndexOffset = viOffset || 0;
	this._voxelsOffset = voOffset || 0;
	this._children = [];
	VoxelGrid.count++;
};

VoxelGrid.count = 0;

VoxelGrid2.prototype.createShortCuts = function () {
	// return;
	const size = this._size;
	this._shortcuts = {};
	for (var z = 0; z < size; z += size - 1) for (var y = 0; y < size; y++) for (var x = 0; x < size; x++) {
		var c = z * size * size + y * size + x;
		this._shortcuts[c] = this._shortcuts[c] || {};
		for (var tz = 0; tz < size; tz += size - 1) for (var ty = 0; ty < size; ty++) for (var tx = 0; tx < size; tx++) {
			var tc = tz * size * size + ty * size + tx;
			this._shortcuts[c][tc] = this.findVoxel(x, y, z, tx, ty, tz);
		}
	}
	for (var z = 0; z < size; z++) for (var y = 0; y < size; y += size - 1) for (var x = 0; x < size; x++) {
		var c = z * size * size + y * size + x;
		this._shortcuts[c] = this._shortcuts[c] || {};
		for (var tz = 0; tz < size; tz++) for (var ty = 0; ty < size; ty += size - 1) for (var tx = 0; tx < size; tx++) {
			var tc = tz * size * size + ty * size + tx;
			this._shortcuts[c][tc] = this.findVoxel(x, y, z, tx, ty, tz);
		}
	}
	for (var z = 0; z < size; z++) for (var y = 0; y < size; y++) for (var x = 0; x < size; x += size - 1) {
		var c = z * size * size + y * size + x;
		this._shortcuts[c] = this._shortcuts[c] || {};
		for (var tz = 0; tz < size; tz++) for (var ty = 0; ty < size; ty++) for (var tx = 0; tx < size; tx += size - 1) {
			var tc = tz * size * size + ty * size + tx;
			this._shortcuts[c][tc] = this.findVoxel(x, y, z, tx, ty, tz);
		}
	}
};

VoxelGrid2.prototype.findVoxel = function (x, y, z, tx, ty, tz) {
	const c = vec3(x,y,z);
	const cf = addS(c, 0.5);
	const d = normalize(vec3(tx-x, ty-y, tz-z));
	const step = sign(d);
	const deltaDist = vec3(
		length(mulS(d, 1 / d.x)),
		length(mulS(d, 1 / d.y)),
		length(mulS(d, 1 / d.z))
	);
	const next = mul(add(maxV(vec3(0), step), add(mulS(mul(step, cf), -1), mul(step, c))), deltaDist);
	const size = this._size;
	var ci = c.z * size * size + c.y * size + c.x;
	while (!(c.x < 0 || c.y < 0 || c.z < 0 || c.x >= size || c.y >= size || c.z >= size)) {
		if ((this._voxelIndex[ci >> 3] >> (ci & 7)) & 1) {
			return c;
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

VoxelGrid2.prototype.add = function (obj) {
	var size = this._size;
	var bbox = obj.getBoundingBox();
	var scaleVec = vec3(this._scale);
	var start = maxV(vec3(0), minV(vec3(this._size - 1), floorV(this.toGrid(bbox.min)))); //sub(bbox.min, mulS(this._dims, (this._levelSizeIndex===0?1:0) * 1.01/size))))));
	var end = maxV(vec3(0), minV(vec3(this._size - 1), floorV(this.toGrid(bbox.max)))); //(add(bbox.max, mulS(this._dims, (this._levelSizeIndex===0?1:0) * 1.01/size))))));
	for (var z = start.z; z <= end.z; z++) {
		var cz = z * size * size;
		for (var y = start.y; y <= end.y; y++) {
			var cy = cz + y * size;
			for (var x = start.x; x <= end.x; x++) {
				var c = cy + x;
				bbox.min = this.fromGrid(vec3(x,y,z));
				bbox.max = add(bbox.min, scaleVec);
				if (obj.intersectBox(bbox)) {
					var vi = (this._voxelIndex[c >> 5] >> (c & 31)) & 1;
					if (!vi) {
						this._voxelIndex[c >> 5] |= (1 << (c & 31));
						this._voxels[c] = this._children.length;
						var v = this._leaf === 0 ? new Voxel() : new VoxelGrid2(this.fromGrid(vec3(x, y, z)), scaleVec, this._levelSizes, this._levelSizeIndex + 1);
						this._children.push(v);
					}
					var v = this._children[this._voxels[c]];
					v.add(obj);
				}
			}
		}
	}
};

VoxelGrid2.pruned = 0;
VoxelGrid2.prototype.prune = function() {
	this._children = this._children.map(c => c.prune ? c.prune() : c);
	if (this._children.every(c => c._objects && c._objects.length === 1 && c._objects[0] === this._children[0]._objects[0])) {
		VoxelGrid2.pruned += this._children.length;
		return this._children[0];	
	}
	return this;
};

VoxelGrid2.prototype.fromGrid = function (coord) {
	return add(this._origin, mulS(coord, this._scale));
};

VoxelGrid2.prototype.toGrid = function (point) {
	return divS(sub(point, this._origin), this._scale);
};

VoxelGrid2.prototype.intersect = function (ray) {
	VoxelGrid.descendCount++;
	const t = this._box.intersect(ray);
	if (t === null) {
		return null;
	}

	const cf = this.toGrid(add(ray.o, mulS(ray.d, t.distance + this._scale * 0.0001)));
	const c = floorV(cf);

	const size = this._size;
	var ci = c.z * size * size + c.y * size + c.x;

	if (this._shortcuts) {
		if (!this._shortcuts[ci]) {
			return null;
		} else {
			const eray = new Ray(tray.o, neg(tray.d), tray.time);
			eray.o = add(eray.o, mulS(eray.d, -2 * size));
			eray.o = add(eray.o, mulS(eray.d, this._box.intersect(eray).distance + this._scale * 0.0001));
			const exit = floorV(this.toGrid(eray.o));
			const ei = exit.z * size * size + exit.y * size + exit.x;
			const nc = this._shortcuts[ci][ei];
			if (!nc) {
				return null;
			}
			cf.x += nc.x - c.x;
			cf.y += nc.y - c.y;
			cf.z += nc.z - c.z;
			c.x = nc.x;
			c.y = nc.y;
			c.z = nc.z;
			ci = c.z * size * size + c.y * size + c.x;
		}
	}

	const deltaDist = vec3();
	const step = sign(ray.d);
	deltaDist.x = length(mulS(ray.d, 1 / ray.d.x));
	deltaDist.y = length(mulS(ray.d, 1 / ray.d.y));
	deltaDist.z = length(mulS(ray.d, 1 / ray.d.z));
	const next = mul(add(maxV(vec3(0), step), add(mulS(mul(step, cf), -1), mul(step, c))), deltaDist);
	while (!(c.x < 0 || c.y < 0 || c.z < 0 || c.x >= size || c.y >= size || c.z >= size)) {
		VoxelGrid.stepCount++;
		const vi = (this._voxelIndex[ci >> 5] >> (ci & 31)) & 1;
		if (vi) {
			ray.light.y += 0.01;
			const voxel = this._children[this._voxels[ci]];
			if (voxel._objects) {
				const nml = voxel._objects[0]._nml;
				return { distance: t.distance + this._scale, obj: {color: () => vec3(0.8), normal: () => nml }};
			}
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
};
