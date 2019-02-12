var TransformNode = function(obj, transform) {
	this.object = obj;
	this.transform = transform || new THREE.Object3D();
};

TransformNode.prototype.intersect = function(ray) {
	var inv = this.transform.inverseMatrix;

	var rayO = new THREE.Vector3().copy(ray.o).applyMatrix4(inv);
	var rayD = new THREE.Vector3().copy(ray.d).applyMatrix3(inv).normalize();
	var hit = this.object.intersect(new Ray(rayO, rayD, ray.time));
	if (hit) {
		hit.obj = new TransformNode(hit.obj, this.transform);
	}
	return hit;
};

TransformNode.prototype.normal = function(point) {
	var inv = this.transform.normalMatrix;

	var tp = new THREE.Vector3().copy(point).applyMatrix4(this.transform.matrix);

	var nml = this.object.normal(tp);
	return new THREE.Vector3().copy(nml).applyMatrix3(inv);
};

TransformNode.prototype.color = function() {
	return this.object.color();
};

TransformNode.prototype.getBoundingBox = function() {
	var inv = this.transform.matrix;

	var box = this.object.getBoundingBox();	
    var boxVertices = [
    	box.min, 
    	box.max,
    	vec3(box.min.x, box.min.y, box.max.z),
    	vec3(box.min.x, box.max.y, box.min.z),
    	vec3(box.min.x, box.max.y, box.max.z),
    	vec3(box.max.x, box.max.y, box.min.z),
    	vec3(box.max.x, box.min.y, box.max.z),
    	vec3(box.max.x, box.min.y, box.min.z)
    ];
    var min = vec3(1/0);
    var max = vec3(-1/0);
    for (var i=0; i<boxVertices.length; i++) {
    	var v = boxVertices[i];
    	var tv = new THREE.Vector3().copy(v).applyMatrix4(inv);
    	min = minV(min, tv);
    	max = maxV(max, tv);
    }
    return new Box(min, max);
};