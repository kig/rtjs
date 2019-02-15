var BVHNode = function(objects, nodeSize, axis) {
	if (nodeSize === undefined) {
		nodeSize = 2;
	}
	if (axis === undefined) {
		axis = 0;
	}
	axis = axis % 3;
	this.axis = axis;
	this.nodeSize = nodeSize;
	this.subTreeSize = 1;
	if (axis === 0) {
		objects.sort(BVHNode.xCompare);
	} else if (axis === 1) {
		objects.sort(BVHNode.yCompare);
	} else {
		objects.sort(BVHNode.zCompare);
	}
	this.objects = objects;
	this.leaf = false;
	if (objects.length <= 2) {
		this.leaf = true;
		this.children = objects;
	} else {
		if (objects.length <= this.nodeSize) {
			this.nodeSize = 2;
		}
		this.children = [];
		var chunkSize = ceil(objects.length / this.nodeSize);
		for (var i=0; i<objects.length; i+=chunkSize) {
			var node = new BVHNode(objects.slice(i, i+chunkSize), this.nodeSize, axis + 1)
			this.children.push(node);
			this.subTreeSize += node.subTreeSize;
		}
	}
	var min = vec3(1/0);
	var max = vec3(-1/0);
	for (var i=0; i<this.children.length; i++) {
		var node = this.children[i];
		min = minV(min, node.getBoundingBox().min),
		max = maxV(max, node.getBoundingBox().max)
	};
	this._box = new Box(min, max);
};


BVHNode.visitedCount = 0;
BVHNode.primitiveTests = 0;

BVHNode.prototype.intersect = function(ray, closestHit, bh) {
	if (closestHit === undefined) {
		closestHit = {obj: null, distance: Infinity};
	}
	BVHNode.visitedCount++;
	var t = bh || this._box.intersect(ray);
	if (t && t.distance < closestHit.distance) {
		const axis = this.axis;
		let backwards = false;
		if (axis === 0 && ray.d.x < 0) {
			backwards = true;
		} else if (axis === 1 && ray.d.y < 0) {
			backwards = true;
		} else if (ray.d.z < 0) {
			backwards = true;
		}
		var leaf = this.leaf ? 1 : 0;
		BVHNode.primitiveTests += this.children.length * leaf;
		if (backwards) {
			for (var i=0; i<this.children.length; i++) {
				if (!this.leaf) {
					var bh = this.children[i]._box.intersect(ray);
					if (!bh) {
						continue;
					}
					if (bh.distance >= closestHit.distance) {
						break;
					}
					var tc = this.children[i].intersect(ray, closestHit, bh);
				} else {
					var tc = this.children[i].intersect(ray, closestHit);
				}
				if (tc && tc.distance < closestHit.distance) {
					closestHit.obj = tc.obj;
					closestHit.distance = tc.distance;
				}
			}
		} else {
			for (var i=this.children.length-1; i>=0; i--) {
				if (!this.leaf) {
					var bh = this.children[i]._box.intersect(ray);
					if (!bh) {
						continue;
					}
					if (bh.distance >= closestHit.distance) {
						break;
					}
					var tc = this.children[i].intersect(ray, closestHit, bh);
				} else {
					var tc = this.children[i].intersect(ray, closestHit);
				}
				if (tc && tc.distance < closestHit.distance) {
					closestHit.obj = tc.obj;
					closestHit.distance = tc.distance;
				}
			}
		}
	}
	return closestHit.obj ? closestHit : null;
};

BVHNode.prototype.getBoundingBox = function() {
	return this._box;
};

BVHNode.xCompare = function(a, b) {
	return a.getBoundingBox().min.x - b.getBoundingBox().min.x;
};

BVHNode.yCompare = function(a, b) {
	return a.getBoundingBox().min.y - b.getBoundingBox().min.y;
};

BVHNode.zCompare = function(a, b) {
	return a.getBoundingBox().min.z - b.getBoundingBox().min.z;
};

var FastBVH = function(objects) {
	var bvh = new Uint16Array(ceil(objects.length/2)*2);
	var boxes = [];
	var axis = 0;
	var start = 0;
	var mid = floor(objects.length/2);
	var end = objects.length;
};
