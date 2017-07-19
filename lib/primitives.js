
var Sphere = function(center, radius, color) {
	this.center = center;
	this.radius = radius;
	this._color = color;
	this._lastIntersect = {o: null, d: null, t: -1/0};
};

Sphere.prototype.intersect = function(ray) {
	if (this._lastIntersect.o === ray.o &&
		this._lastIntersect.d === ray.d) {
		return this._lastIntersect.t;
	}
	this._lastIntersect.o = ray.o;
	this._lastIntersect.d = ray.d;
	var rc = sub(ray.o, this.center); 
	var c = dot(rc, rc) - this.radius*this.radius;
	var b = dot(ray.d, rc);
	var d = b*b - c;
	if (d < 0) {
		this._lastIntersect.t = -1/0;
		return -1/0;
	}
	var t = -b - sqrt(d);
	this._lastIntersect.t = t;
	return t;
};

Sphere.prototype.normal = function(point) {
	return normalize(sub(point, this.center));
};

Sphere.prototype.color = function(ray) {
	return this._color;
};

Sphere.prototype.getBoundingBox = function() {
	return new Box(
		sub(this.center, vec3(this.radius)),
		add(this.center, vec3(this.radius))
	);
};

Sphere.prototype.intersectBox = function(box) {
	var check = function(v, min, max) {
		if (v < min) {
			var val = min - v;
			return val * val;
		}
		if (v > max) {
			var val = v - max;
			return val * val;
		}
		return 0;
	};
	var sq = 0;
	var c = this.center;
	sq += check(c.x, box.min.x, box.max.x);
	sq += check(c.y, box.min.y, box.max.y);
	sq += check(c.z, box.min.z, box.max.z);
	return sq <= this.radius * this.radius;
};

var Plane = function(point, normal, color) {
	this._point = point;
	this._normal = normalize(normal);
	this._color = color;
};

Plane.prototype.intersect = function(ray) {
	var pd = dot(this._normal, ray.d);
	if (abs(pd) > 0.00001) {
		var dist = dot(this._normal, sub(this._point, ray.o)) / pd;
		return dist;
	}
	return -1/0;
};

Plane.prototype.intersectInf = function(ray) {
	var pd = dot(this._normal, ray.d);
	if (abs(pd) > 0.00001) {
		var dist = dot(this._normal, sub(this._point, ray.o)) / pd;
		return dist;
	}
	return 1/0;
};

Plane.prototype.normal = function(point) {
	return this._normal;
};

Plane.prototype.color = function(ray) {
	if (abs((floor(ray.o.x*2) + floor(ray.o.z*2)) % 2) == 0) {
		return vec3(0.2);
	}
	return this._color;
};

Plane.prototype.getBoundingBox = function() {
	return {
		min: vec3(-1/0, -1/0, -1/0),
		max: vec3(1/0, 1/0, 1/0)
	};
};

var Triangle = function(vertices, normals, color) {
	this._vertices = vertices;
	this._normals = normals;
	this._color = color;
	this._lastIntersect = {o: null, d: null, t: -1/0};
    this._nml = cross(
    	sub(this._vertices[1], this._vertices[0]), 
    	sub(this._vertices[2], this._vertices[0])
	);
    this._edges = [
        sub(this._vertices[0], this._vertices[1]),
        sub(this._vertices[1], this._vertices[2]),
        sub(this._vertices[2], this._vertices[0]),
    ];
    this._centroid = mulS(add(add(vertices[0], vertices[1]), vertices[2]), 1/3);
};

Triangle.prototype.intersect = function(ray) {
	if (this._lastIntersect.o === ray.o &&
		this._lastIntersect.d === ray.d) {
		return this._lastIntersect.t;
	}
	this._lastIntersect.o = ray.o;
	this._lastIntersect.d = ray.d;

	var e1,e2,h,s,q; // vec3
	var a,f,u,v,t; // float

	var v0 = this._vertices[0];
	var v1 = this._vertices[1];
	var v2 = this._vertices[2];

	e1 = sub(v1,v0);
	e2 = sub(v2,v0);

	h = cross(ray.d,e2);
	a = dot(e1,h);

	f = 1.0 / a;
	s = sub(ray.o, v0);
	u = f * dot(s,h);

	q = cross(s,e1);
	v = f * dot(ray.d,q);

	if (u < 0.0 || u > 1.0 || v < 0.0 || u+v > 1.0) {
		this._lastIntersect.t = -1/0;
		return -1/0;
	}

	// at this stage we can compute t to find out where
	// the intersection point is on the line
	t = f * dot(e2,q);

	if (t <= 0.00001 || (a > -0.00001 && a < 0.00001)) {
		this._lastIntersect.t = -1/0;
		return -1/0;
	}
	this._lastIntersect.t = t;
	return t;
};

Triangle.prototype.normal = function(point) {
	var v0 = this._vertices[0];
	var v1 = this._vertices[1];
	var v2 = this._vertices[2];
	var u = dot(point, sub(v1, v0));
	var v = dot(point, sub(v2, v0));
	var n0 = this._normals[0];
	var n1 = this._normals[1];
	var n2 = this._normals[2];
	return normalize(mulS(add(mix(n0, n1, u), mix(n0, n2, v)), 0.5));
};

Triangle.prototype.color = function(ray) {
	return this._color;
};

Triangle.prototype.getBoundingBox = function() {
	return new Box(
		vec3(
			min(this._vertices[0].x, this._vertices[1].x, this._vertices[2].x),
			min(this._vertices[0].y, this._vertices[1].y, this._vertices[2].y),
			min(this._vertices[0].z, this._vertices[1].z, this._vertices[2].z)
			),
		vec3(
			max(this._vertices[0].x, this._vertices[1].x, this._vertices[2].x),
			max(this._vertices[0].y, this._vertices[1].y, this._vertices[2].y),
			max(this._vertices[0].z, this._vertices[1].z, this._vertices[2].z)
			)
	);
};

Triangle.prototype.intersectBoxProject = function(points, axis, minMax) {
    var min = 1/0;
    var max = -1/0;
    for (var i=0; i<points.length; i++) {
        var val = dot(axis, points[i]);
        if (val < min) min = val;
        if (val > max) max = val;
    }
    minMax.min = min;
    minMax.max = max;
};

Triangle.boxNormals = [vec3(1,0,0), vec3(0,1,0), vec3(0,0,1)];

Triangle.prototype.intersectBox = function(box) {
	if (box.contains(this._centroid) || box.contains(this._vertices[0]) || box.contains(this._vertices[1]) || box.contains(this._vertices[2])) {
		return true;
	}

    var boxNormals = Triangle.boxNormals;
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

    // Test the box normals (x-, y- and z-axes)
    var minMax = {min: 0, max: 0};
    var triMinMax = {min: 0, max: 0};
    this.intersectBoxProject(this._vertices, boxNormals[0], minMax);
    if (minMax.max < box.min.x || minMax.min > box.max.x)
        return false; // No intersection possible.
    this.intersectBoxProject(this._vertices, boxNormals[1], minMax);
    if (minMax.max < box.min.y || minMax.min > box.max.y)
        return false; // No intersection possible.
    this.intersectBoxProject(this._vertices, boxNormals[2], minMax);
    if (minMax.max < box.min.z || minMax.min > box.max.z)
        return false; // No intersection possible.

    // Test the triangle normal
    var triangleOffset = dot(this._nml, this._vertices[0]);
    this.intersectBoxProject(boxVertices, this._nml, minMax);
    if (minMax.max < triangleOffset || minMax.min > triangleOffset)
        return false; // No intersection possible.

    // Test the nine edge cross-products
    var triangleEdges = this._edges;
    for (var i = 0; i < 3; i++)
    for (var j = 0; j < 3; j++)
    {
        // The box normals are the same as it's edge tangents
        var axis = cross(triangleEdges[i], boxNormals[j]);
        this.intersectBoxProject(boxVertices, axis, minMax);
        this.intersectBoxProject(this._vertices, axis, triMinMax);
        if (minMax.max < triMinMax.min || minMax.min > triMinMax.max)
            return false; // No intersection possible
    }

    // No separating axis found.
	return true;
};

var Box = function(min, max) {
	this.min = min;
	this.max = max;
};

Box.prototype.contains = function(p) {
	return (
		this.min.x <= p.x && this.min.y <= p.y && this.min.z <= p.z &&
		this.max.x >= p.x && this.max.y >= p.y && this.max.z >= p.z
	);
};

Box.prototype.intersect = function(r) {
	var tmin = (this.min.x - r.o.x) / r.d.x; 
	var tmax = (this.max.x - r.o.x) / r.d.x; 

	if (tmin > tmax) {
		var tmp = tmin;
		tmin = tmax;
		tmax = tmp;
	}

	var tymin = (this.min.y - r.o.y) / r.d.y; 
	var tymax = (this.max.y - r.o.y) / r.d.y; 

	if (tymin > tymax) {
		var tmp = tymin;
		tymin = tymax;
		tymax = tmp;
	}

	if ((tmin > tymax) || (tymin > tmax)) 
		return -1/0; 

	if (tymin > tmin) 
		tmin = tymin;

	if (tymax < tmax) 
		tmax = tymax; 

	var tzmin = (this.min.z - r.o.z) / r.d.z; 
	var tzmax = (this.max.z - r.o.z) / r.d.z; 

	if (tzmin > tzmax) {
		var tmp = tzmin;
		tzmin = tzmax;
		tzmax = tmp;
	}

	if ((tmin > tzmax) || (tzmin > tmax)) 
		return -1/0; 

	if (tzmin > tmin) 
		tmin = tzmin; 

	if (tzmax < tmax) 
		tmax = tzmax; 

	return tmin;
};
