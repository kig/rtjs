
var Voxel = function() {
	this._objects = [];
};

Voxel.prototype.add = function(obj) {
	this._objects.push(obj);
};

Voxel.prototype.intersect = function(ray) {
	VoxelGrid.cmpCount += this._objects.length;
	return intersect(ray, this._objects);
};


var VoxelGrid = function(size, origin, dims, leaf, levelSize) {
	this._origin = origin;
	this._dims = dims;
	this._size = size;
	this._scale = this._dims.x / this._size;
	this._leaf = leaf;
	this._levelSize = levelSize;
	var voxels = [];
	for (var z=0; z<size; z++) {
		var vz = voxels[z] = [];
		for (var y=0; y<size; y++) {
			var vy = vz[y] = [];
			for (var x=0; x<size; x++) {
				vy[x] = null;
			}
		}
	}
	this._voxels = voxels;
};
VoxelGrid.stepCount = 0;
VoxelGrid.cmpCount = 0;

VoxelGrid.prototype.add = function(obj) {
	var bbox = obj.getBoundingBox();
	var start = maxV(vec3(0), minV(vec3(this._size-1), floorV(this.toGrid(bbox.min))));
	var end = maxV(vec3(0), minV(vec3(this._size-1), floorV(this.toGrid(bbox.max))));
	var box = new Box;
	for (var z=start.z; z<=end.z; z++) {
		for (var y=start.y; y<=end.y; y++) {
			for (var x=start.x; x<=end.x; x++) {
				box.min = this.fromGrid(vec3(x, y, z));
				box.max = this.fromGrid(vec3(x+1, y+1, z+1));
				if (obj.intersectBox(box)) {
					var v = this._voxels[z][y][x];
					if (!v) {
						v = this._voxels[z][y][x] = this._leaf === 0 ? new Voxel() : new VoxelGrid(this._levelSize, this.fromGrid(vec3(x,y,z)), vec3(this._scale), this._leaf-1, this._levelSize);
					}
					v.add(obj);
				}
			}
		}
	}
};

VoxelGrid.prototype.fromGrid = function(coord) {
	return add(this._origin, mulS(coord, this._scale));
};

VoxelGrid.prototype.toGrid = function(point) {
	return mulS(sub(point, this._origin), 1/this._scale);
};

VoxelGrid.prototype.intersect = function(ray) {
	var tray = new Ray(ray.o, ray.d, ray.time);
	//debugger;
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
	/*
	var next = mul(vec3(
		step.x < 0 ? (cf.x - c.x) : (c.x + 1.0 - cf.x),
		step.y < 0 ? (cf.y - c.y) : (c.y + 1.0 - cf.y),
		step.z < 0 ? (cf.z - c.z) : (c.z + 1.0 - cf.z)
	), deltaDist);
	*/
	while (!(any(lessThan(c, vec3(0))) || any(greaterThanEqual(c, vec3(this._size))))) {
		VoxelGrid.stepCount++;
		var voxel = this._voxels[c.z][c.y][c.x];
		if (voxel) {
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
			} else {
				next.z += deltaDist.z;
				c.z += step.z;
			}
		} else if (next.y < next.z) {
			next.y += deltaDist.y;
			c.y += step.y;
		} else {
			next.z += deltaDist.z;
			c.z += step.z;
		}
	}
	return null;
};
