var VoxelGrid2 = function(size, origin, dims, leaf, levelSize) {
	this._origin = origin;
	this._dims = dims;
	this._size = size;
	this._scale = this._dims.x / this._size;
	this._leaf = leaf;
	this._levelSize = levelSize;
	this._box = new Box(origin, add(origin, dims));
	this._voxels = new Uint16Array(size*size*size);
	this._children = [];
};

VoxelGrid2.prototype.add = function(obj) {
	var size = this._size;
	var bbox = obj.getBoundingBox();
	var start = maxV(vec3(0), minV(vec3(this._size-1), floorV(this.toGrid(bbox.min))));
	var end = maxV(vec3(0), minV(vec3(this._size-1), floorV(this.toGrid(bbox.max))));
	var box = new Box;
	for (var z=start.z; z<=end.z; z++) {
		var cz = z * size * size;
		for (var y=start.y; y<=end.y; y++) {
			var cy = cz + y*size;
			for (var x=start.x; x<=end.x; x++) {
				box.min = this.fromGrid(vec3(x, y, z));
				box.max = this.fromGrid(vec3(x+1, y+1, z+1));
				if (obj.intersectBox(box)) {
					var c = cy + x;
					var vi = this._voxels[c];
					if (!vi) {
						vi = this._voxels[c] = this._children.length+1;
						var v = this._leaf === 0 ? new Voxel() : new VoxelGrid2(this._levelSize, this.fromGrid(vec3(x,y,z)), vec3(this._scale), this._leaf-1, this._levelSize);
						this._children.push(v);
					}
					var v = this._children[vi-1];
					v.add(obj);
				}
			}
		}
	}
};

VoxelGrid2.prototype.fromGrid = function(coord) {
	return add(this._origin, mulS(coord, this._scale));
};

VoxelGrid2.prototype.toGrid = function(point) {
	return divS(sub(point, this._origin), this._scale);
};

VoxelGrid2.prototype.intersect = function(ray) {
	var tray = new Ray(ray.o, ray.d, ray.time);
	var t = new Box(this._origin, add(this._origin, this._dims)).intersect(tray);

	tray.o = add(tray.o, mulS(tray.d, t+this._scale*0.0001));

	var cf = this.toGrid(tray.o);
	var c = floorV(cf);

	var deltaDist = vec3();
	var step = sign(ray.d);
	deltaDist.x = length(mulS(ray.d, 1/ray.d.x));
	deltaDist.y = length(mulS(ray.d, 1/ray.d.y));
	deltaDist.z = length(mulS(ray.d, 1/ray.d.z));
	var next = mul(add(maxV(vec3(0), step), add(mulS(mul(step, cf), -1), mul(step, c))), deltaDist);
	var size = this._size;
	var ci = c.z * size * size + c.y * size + c.x;
	while (!(c.x < 0 || c.y < 0 || c.z < 0 || c.x >= size || c.y >= size || c.z >= size)) {
		VoxelGrid.stepCount++;
		var vi = this._voxels[ci];
		if (vi) {
			var voxel = this._children[vi-1];
			var hit = voxel.intersect(ray);
			if (hit) {
				var p = add(ray.o, mulS(ray.d, hit.distance));
				if (all(equal(c, floorV(this.toGrid(p))))) {
					return hit;
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
