var BVHNode = function(objects) {
	this.subTreeSize = 1;
	var axis = floor(random() * 3);
	if (axis === 0) {
		objects.sort(BVHNode.xCompare);
	} else if (axis === 1) {
		objects.sort(BVHNode.yCompare);
	} else {
		objects.sort(BVHNode.zCompare);
	}
	this.leaf = false;
	if (objects.length === 1) {
		this.leaf = true;
		this.left = this.right = objects[0];
	} else if (objects.length === 2) {
		this.leaf = true;
		this.left = objects[0];
		this.right = objects[1];
	} else {
		this.left = new BVHNode(objects.slice(0, objects.length/2));
		this.right = new BVHNode(objects.slice(objects.length/2));
		this.subTreeSize += this.left.subTreeSize + this.right.subTreeSize;
	}
	this._box = new Box(
		minV(this.left.getBoundingBox().min, this.right.getBoundingBox().min),
		maxV(this.left.getBoundingBox().max, this.right.getBoundingBox().max)
	);
};

BVHNode.visitedCount = 0;
BVHNode.primitiveTests = 0;

BVHNode.prototype.intersect = function(ray) {
	BVHNode.visitedCount++;
	var t = this.getBoundingBox().intersect(ray);
	if (t > 0) {
		if (this.leaf) {
			var t_left = this.left.intersect(ray);
			var t_right = this.right.intersect(ray);
			if (this.left !== this.right) {
				BVHNode.primitiveTests++;
			}
			BVHNode.primitiveTests++;
			if (t_left > 0 && t_right > 0) {
				if (t_left < t_right) {
					return {obj: this.left, distance: t_left};
				} else {
					return {obj: this.right, distance: t_right};
				}
			} else if (t_left > 0) {
				return {obj: this.left, distance: t_left};
			} else if (t_right > 0) {
				return {obj: this.right, distance: t_right};
			}
		} else {
			var t_left = this.left.intersect(ray);
			var t_right = this.right.intersect(ray);
			if (t_left && t_right) {
				if (t_left.distance < t_right.distance) {
					return t_left;
				} else {
					return t_right;
				}
			} else {
				return t_left || t_right;
			}
		}
	}
	return null;
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
