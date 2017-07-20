
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
	this._box = new Box(origin, add(origin, dims));
	var voxels = [];
	var distances = [];
	for (var z=0; z<size; z++) {
		var vz = voxels[z] = [];
		var dz = distances[z] = [];
		for (var y=0; y<size; y++) {
			var vy = vz[y] = [];
			var dy = dz[y] = [];
			for (var x=0; x<size; x++) {
				vy[x] = null;
				dy[x] = this._size;
			}
		}
	}
	this._voxels = voxels;
	this._distances = distances;
};
VoxelGrid.stepCount = 0;
VoxelGrid.cmpCount = 0;

VoxelGrid.prototype.computeDistanceStep = function() {
	var voxels = this._voxels;
	var distances = this._distances;
	var size = this._size;
	for (var z=0; z<size; z++) {
		var vz = voxels[z];
		var dz = distances[z];
		for (var y=0; y<size; y++) {
			var vy = vz[y];
			var dy = dz[y];
			for (var x=0; x<size; x++) {
				if (vy[x]) {
					dy[x] = 0;
				} else {
					var d = dy[x];
					for (var zz=-1; zz<=1; zz++) {
						if (z + zz < 0 || z + zz >= this._size) {
							continue;
						}
						for (var yy=-1; yy<=1; yy++) {
							if (y + yy < 0 || y + yy >= this._size) {
								continue;
							}
							for (var xx=-1; xx<=1; xx++) {
								if (x + xx < 0 || x + xx >= this._size) {
									continue;
								}
								if (xx === yy && yy === zz && zz === 0) {
									continue;
								}
								d = min(d, distances[z+zz][y+yy][x+xx] + 1);
							}
						}
					}
					dy[x] = d;
				}
			}
		}
	}
};

VoxelGrid.prototype.computeDistances = function() {
	for (var i=0; i<this._size; i++) {
		this.computeDistanceStep();
	}
	if (this._leaf > 0) {
		var voxels = this._voxels;
		var size = this._size;
		for (var z=0; z<size; z++) {
			var vz = voxels[z];
			for (var y=0; y<size; y++) {
				var vy = vz[y];
				for (var x=0; x<size; x++) {
					if (vy[x]) {
						vy[x].computeDistances();
					}
				}
			}
		}
	}
};

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
	return divS(sub(point, this._origin), this._scale);
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
	var size = this._size;
	while (!(c.x < 0 || c.y < 0 || c.z < 0 || c.x >= size || c.y >= size || c.z >= size)) {
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