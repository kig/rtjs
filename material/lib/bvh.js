/*
	BVH splits the scene into possibly overlapping bounding boxes.
	Because of the potential overlaps, BVH traversal doesn't terminate on first found primitive.
*/

var BVHNode = function(objects, nodeSize, axis) {
	if (nodeSize === undefined) {
		nodeSize = 2;
	}
	if (axis === undefined) {
		axis = 0;
		objects.forEach(element => {
			element._box = null;
		});
	}
	axis = axis % 3;
	this.axis = axis;
	this.nodeSize = nodeSize;
	this.subTreeSize = 1;
	this.totalArea = objects.reduce((s, o) => s + o.area, 0);
	if (axis === 0) {
		objects.sort(BVHNode.xCompare);
	} else if (axis === 1) {
		objects.sort(BVHNode.yCompare);
	} else {
		objects.sort(BVHNode.zCompare);
	}
	this.leaf = false;
	if (objects.length < 2) {
		this.leaf = true;
		this.children = objects;
	} else {
		if (objects.length <= this.nodeSize) {
			this.nodeSize = 2;
		}
		this.children = [];
		var useSAH = true;
		if (useSAH) {
			var chunkAreaTarget = this.totalArea / this.nodeSize;
			var chunkArea = 0;
			var chunkStartIdx = 0;
			for (var i = 0; i < objects.length-1; i++) {
				chunkArea += objects[i].area;
				if (chunkArea >= chunkAreaTarget) {
					var node = new BVHNode(objects.slice(chunkStartIdx, i+1), this.nodeSize, axis + 1)
					this.children.push(node);
					this.subTreeSize += node.subTreeSize;
					chunkStartIdx = i+1;
					chunkArea = 0;
				}
			}
			if (chunkStartIdx === 0) {
				var node = new BVHNode(objects.slice(chunkStartIdx, objects.length-1), this.nodeSize, axis + 1)
				this.children.push(node);
				this.subTreeSize += node.subTreeSize;
				var node = new BVHNode(objects.slice(objects.length-1), this.nodeSize, axis + 1)
				this.children.push(node);
				this.subTreeSize += node.subTreeSize;
			} else {
				var node = new BVHNode(objects.slice(chunkStartIdx), this.nodeSize, axis + 1)
				this.children.push(node);
				this.subTreeSize += node.subTreeSize;
			}
		} else {
			var chunkSize = ceil(objects.length / this.nodeSize);
			for (var i=0; i<objects.length; i+=chunkSize) {
				var node = new BVHNode(objects.slice(i, i+chunkSize), this.nodeSize, axis + 1)
				this.children.push(node);
				this.subTreeSize += node.subTreeSize;
			}
		}
	}
	var bmin = vec3(1/0);
	var bmax = vec3(-1/0);
	for (var i=0; i<this.children.length; i++) {
		var node = this.children[i];
		var nbox = node.getBoundingBox();
		bmin = minV(bmin, nbox.min);
		bmax = maxV(bmax, nbox.max);
	};
	this._box = new Box(bmin, bmax);
};

BVHNode.visitedCount = 0;
BVHNode.primitiveTests = 0;

BVHNode.prototype.intersect = function(ray, closestHit, boxHit) {
	if (closestHit === undefined) {
		closestHit = {obj: null, distance: 1e9};
	}
	BVHNode.visitedCount++;
	const axis = this.axis;
	let backwards = false;
	if (axis === 0 && ray.d.x < 0) {
		backwards = true;
	} else if (axis === 1 && ray.d.y < 0) {
		backwards = true;
	} else if (axis === 2 && ray.d.z < 0) {
		backwards = true;
	}
	if (this.leaf) {
		BVHNode.primitiveTests += this.children.length;
		ray.light.x += 0.1;
	} else {
		ray.light.z += 0.1;
	}
	if (!backwards) {
		for (var i=0; i<this.children.length; i++) {
			var boxHit = this.children[i]._box.intersect(ray);
			if (!boxHit) continue;
			if (boxHit.distance >= closestHit.distance) break;
			var tc = this.children[i].intersect(ray, closestHit);
			if (tc && tc.distance < closestHit.distance) {
				closestHit.obj = tc.obj;
				closestHit.distance = tc.distance;
			}
		}
	} else {
		for (var i=this.children.length-1; i >= 0; i--) {
			var boxHit = this.children[i]._box.intersect(ray);
			if (!boxHit) continue;
			if (boxHit.distance >= closestHit.distance) break;
			var tc = this.children[i].intersect(ray, closestHit);
			if (tc && tc.distance < closestHit.distance) {
				closestHit.obj = tc.obj;
				closestHit.distance = tc.distance;
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
