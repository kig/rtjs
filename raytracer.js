var Vec3 = function(x, y, z) {
	this.x = x;
	this.y = y;
	this.z = z;
};

Object.defineProperties(Vec3.prototype, {
	0: { get: function() { return this.x; }, set: function(v) { this.x = v; }},
	1: { get: function() { return this.y; }, set: function(v) { this.y = v; }},
	2: { get: function() { return this.z; }, set: function(v) { this.z = v; }},

	xxx: { get: function() { return new Vec3(this.x, this.x, this.x); }},
	xxy: { get: function() { return new Vec3(this.x, this.x, this.y); }},
	xxz: { get: function() { return new Vec3(this.x, this.x, this.z); }},
	xyx: { get: function() { return new Vec3(this.x, this.y, this.x); }},
	xyy: { get: function() { return new Vec3(this.x, this.y, this.y); }},
	xyz: { get: function() { return new Vec3(this.x, this.y, this.z); }},
	xzx: { get: function() { return new Vec3(this.x, this.z, this.x); }},
	xzy: { get: function() { return new Vec3(this.x, this.z, this.y); }},
	xzz: { get: function() { return new Vec3(this.x, this.z, this.z); }},

	yxx: { get: function() { return new Vec3(this.y, this.x, this.x); }},
	yxy: { get: function() { return new Vec3(this.y, this.x, this.y); }},
	yxz: { get: function() { return new Vec3(this.y, this.x, this.z); }},
	yyx: { get: function() { return new Vec3(this.y, this.y, this.x); }},
	yyy: { get: function() { return new Vec3(this.y, this.y, this.y); }},
	yyz: { get: function() { return new Vec3(this.y, this.y, this.z); }},
	yzx: { get: function() { return new Vec3(this.y, this.z, this.x); }},
	yzy: { get: function() { return new Vec3(this.y, this.z, this.y); }},
	yzz: { get: function() { return new Vec3(this.y, this.z, this.z); }},

	zxx: { get: function() { return new Vec3(this.z, this.x, this.x); }},
	zxy: { get: function() { return new Vec3(this.z, this.x, this.y); }},
	zxz: { get: function() { return new Vec3(this.z, this.x, this.z); }},
	zyx: { get: function() { return new Vec3(this.z, this.y, this.x); }},
	zyy: { get: function() { return new Vec3(this.z, this.y, this.y); }},
	zyz: { get: function() { return new Vec3(this.z, this.y, this.z); }},
	zzx: { get: function() { return new Vec3(this.z, this.z, this.x); }},
	zzy: { get: function() { return new Vec3(this.z, this.z, this.y); }},
	zzz: { get: function() { return new Vec3(this.z, this.z, this.z); }},
});

var vec3 = function(x, y, z) {
	x = x === undefined ? 0 : x;
	y = y === undefined ? x : y;
	z = z === undefined ? y : z;
	return new Vec3(x,y,z);
};

var sqrt = Math.sqrt;
var abs = Math.abs;
var min = Math.min;
var max = Math.max;
var random = Math.random;
var pow = Math.pow;
var exp = Math.exp;
var sin = Math.sin;
var cos = Math.cos;
var tan = Math.tan;
var floor = Math.floor;
var ceil = Math.ceil;
var round = Math.round;

var sat = function(v) {
	return max(0, min(1, v));
};

var mul = function(u,v) {
	return vec3(u.x*v.x, u.y*v.y, u.z*v.z);
};

var div = function(u,v) {
	return vec3(u.x/v.x, u.y/v.y, u.z/v.z);
};

var add = function(u,v) {
	return vec3(u.x+v.x, u.y+v.y, u.z+v.z);
};

var sub = function(u,v) {
	return vec3(u.x-v.x, u.y-v.y, u.z-v.z);
};

var recip = function(u) {
	return Sdiv(1, u);
};

var mulS = function(u,s) {
	return vec3(u.x*s, u.y*s, u.z*s);
};

var divS = function(u,s) {
	return vec3(u.x/s, u.y/s, u.z/s);
};

var Sdiv = function(s,u) {
	return vec3(s/u.x, s/u.y, s/u.z);
};

var addS = function(u,s) {
	return vec3(u.x+s, u.y+s, u.z+s);
};

var subS = function(u,s) {
	return vec3(u.x-s, u.y-s, u.z-s);
};

var Ssub = function(s,u) {
	return vec3(s-u.x, s-u.y, s-u.z);
};

var dot = function(u,v) {
	return (u.x*v.x + u.y*v.y + u.z*v.z);
};

var reflect = function(v, nml) {
	return sub(v, mulS(nml, 2*dot(v,nml)));
};

var refract = function(v, nml, eta) {
	var d = dot(nml, v);
	var k = 1.0 - eta * eta * (1.0 - d * d);
	if (k < 0.0) {
		return vec3(0);
	}
    return sub(mulS(v, eta), mulS(nml, (eta * d + sqrt(k))));
};

var length = function(v) {
	return sqrt(dot(v,v));
};

var lengthSq = function(v) {
	return dot(v,v);
};

var normalize = function(v) {
	return mulS(v, 1 / length(v));
};

var cross = function(u,v) {
	return vec3(u.y*v.z-u.z*v.y, u.z*v.x-u.x*v.z, u.x*v.y-u.y*v.x);
};

var mix = function(u, v, t) {
	return add(mulS(u, 1-t), mulS(v, t));
};

var neg = function(v) {
	return vec3(-v.x, -v.y, -v.z);
};

var expV = function(v) {
	return vec3(exp(v.x), exp(v.y), exp(v.z));
};

var powV = function(v, u) {
	return vec3(pow(v.x, u.x), pow(v.y, u.y), pow(v.z, u.z));
};

var sqrtV = function(v) {
	return vec3(sqrt(v.x), sqrt(v.y), sqrt(v.z));
};

var floorV = function(v) {
	return vec3(floor(v.x), floor(v.y), floor(v.z));
};

var ceilV = function(v) {
	return vec3(ceil(v.x), ceil(v.y), ceil(v.z));
};

var roundV = function(v) {
	return vec3(round(v.x), round(v.y), round(v.z));
};

var maxV = function(u, v) {
	return vec3(max(u.x, v.x), max(u.y, v.y), max(u.z, v.z));
};

var minV = function(u, v) {
	return vec3(min(u.x, v.x), min(u.y, v.y), min(u.z, v.z));
};

var sign = function(v) {
	return vec3(
		v.x < 0 ? -1 : 1,
		v.y < 0 ? -1 : 1,
		v.z < 0 ? -1 : 1
	);
};

var absV = function(v) {
	return vec3(abs(v.x), abs(v.y), abs(v.z));
};

var any = function(v) {
	return v.x || v.y || v.z;
};

var all = function(v) {
	return v.x && v.y && v.z;
};

var equal = function(u, v) {
	return vec3(u.x === v.x, u.y === v.y, u.z === v.z);
};

var lessThan = function(v, cmp) {
	return vec3(v.x < cmp.x, v.y < cmp.y, v.z < cmp.z);
};

var lessThanEqual = function(v, cmp) {
	return vec3(v.x <= cmp.x, v.y <= cmp.y, v.z <= cmp.z);
};

var greaterThan = function(v, cmp) {
	return vec3(v.x > cmp.x, v.y > cmp.y, v.z > cmp.z);
};

var greaterThanEqual = function(v, cmp) {
	return vec3(v.x >= cmp.x, v.y >= cmp.y, v.z >= cmp.z);
};

var randomVec3Positive = function() {
	return vec3(random(), random(), random());
};

var randomVec3 = function() {
	return subS(mulS(randomVec3Positive(), 2), 1);
};

var randomVec3Unit = function() {
	return normalize(randomVec3());
};

var diskPoint = function() {
	var a = random();
	var r = random();
	return mulS(vec3(cos(a), sin(a), 0.0), sqrt(r));
};

var Ray = function(o, d, time) {
	this.o = o;
	this.d = d;
	this.transmit = vec3(1.0);
	this.light = vec3(0.0);
	this.time = time;
	this.bounce = 0;
};

var intersect = function(ray, scene) {
	var minT = 1/0;
	var obj = null;
	for (var i=0; i<scene.length; i++) {
		var t = scene[i].intersect(ray);
		if (t > 0 && t < minT) {
			minT = t;
			obj = scene[i];
		}
	}
	if (!obj) {
		return null;
	}
	return {distance: minT, obj: obj};
};

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

Triangle.prototype.intersectBox = function(box) {
	var project = function(points, axis, minMax) {
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

    var boxNormals = [vec3(1,0,0), vec3(0,1,0), vec3(0,0,1)];
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
    project(this._vertices, boxNormals[0], minMax);
    if (minMax.max < box.min.x || minMax.min > box.max.x)
        return false; // No intersection possible.
    project(this._vertices, boxNormals[1], minMax);
    if (minMax.max < box.min.y || minMax.min > box.max.y)
        return false; // No intersection possible.
    project(this._vertices, boxNormals[2], minMax);
    if (minMax.max < box.min.z || minMax.min > box.max.z)
        return false; // No intersection possible.

    // Test the triangle normal
    var triangleOffset = dot(this._normals[0], this._vertices[0]);
    project(boxVertices, this._normals[0], minMax);
    if (minMax.max < triangleOffset || minMax.min > triangleOffset)
        return false; // No intersection possible.

    // Test the nine edge cross-products
    var triangleEdges = [
        sub(this._vertices[0], this._vertices[1]),
        sub(this._vertices[1], this._vertices[2]),
        sub(this._vertices[2], this._vertices[0]),
    ];
    for (var i = 0; i < 3; i++)
    for (var j = 0; j < 3; j++)
    {
        // The box normals are the same as it's edge tangents
        var axis = cross(triangleEdges[i], boxNormals[j]);
        project(boxVertices, axis, minMax);
        project(this._vertices, axis, triMinMax);
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

var setupRay = function(camera, x, y, w, h) {
	var uv = new THREE.Vector2(x/w*2 - 1, y/h*2 - 1);

	var origin = new THREE.Vector3();
	origin.setFromMatrixPosition( camera.matrixWorld );

	var direction = new THREE.Vector3();
	direction.set( uv.x, uv.y, 0.5 ).unproject( camera ).sub( origin ).normalize();

	var target = new THREE.Vector3();
	target.addVectors( origin, direction.multiplyScalar(camera.position.distanceTo(camera.target)) );

	origin.add( new THREE.Vector3().copy(diskPoint()).multiplyScalar(0.3).applyMatrix4(camera.matrixWorld) );
	direction.subVectors( target, origin ).normalize();

	var ray = new Ray(vec3(origin.x, origin.y, origin.z), vec3(direction.x, direction.y, direction.z), 0);
	return ray;
};


var loader = new THREE.OBJLoader;
loader.load('bunny.obj', function(bunny) {
	var camera = new THREE.PerspectiveCamera(35, 1, 0.1, 100);
	camera.position.set(3.5, 1.5, 3.5).normalize().multiplyScalar(4.0);
	camera.target = new THREE.Vector3(0, 0.5, 0);
	camera.lookAt(camera.target);
	camera.updateProjectionMatrix();
	camera.updateMatrixWorld();

	var canvasSize = 100;
	var AA_SIZE = 1;

	var canvas = document.createElement('canvas');
	canvas.width = canvas.height = canvasSize;

	document.body.appendChild(canvas);

	var ctx = canvas.getContext('2d');
	var id = ctx.getImageData(0, 0, canvasSize, canvasSize);

	var scene = [
		new Sphere(vec3(2.5,1,1), 1, vec3(1.0, 0.7, 0.3))
	];

	bunny.children[0].geometry.computeBoundingBox();
	var bbox = bunny.children[0].geometry.boundingBox;
	var verts = bunny.children[0].geometry.attributes.position.array;
	var normals = bunny.children[0].geometry.attributes.normal.array;

	var scale = 2.0 / (bbox.max.y - bbox.min.y);
	var xOffset = -(bbox.max.x+bbox.min.x)/2;
	var zOffset = -(bbox.max.z+bbox.min.z)/2;
	var yOffset = -bbox.min.y;

	for (var i = 0; i < verts.length; i += 3) {
		verts[i] += xOffset;
		verts[i] *= scale;
		verts[i+1] += yOffset;
		verts[i+1] *= scale;
		verts[i+2] += zOffset;
		verts[i+2] *= scale;
	}

	var bunnyTris = [];

	for (var i = 0; i < verts.length; i += 3*3) {
		var u = vec3(verts[i], verts[i+1], verts[i+2]);
		var v = vec3(verts[i+3], verts[i+4], verts[i+5]);
		var w = vec3(verts[i+6], verts[i+7], verts[i+8]);
		var x = vec3(normals[i], normals[i+1], normals[i+2]);
		var y = vec3(normals[i+3], normals[i+4], normals[i+5]);
		var z = vec3(normals[i+6], normals[i+7], normals[i+8]);
		bunnyTris.push(new Triangle([u,v,w], [x,y,z], add(vec3(0.8), mulS(absV(u.xxy),0.2))));
	}

	for (var i=0; i<100; i++) {
		var c = mulS(randomVec3Unit(), 4 + random() * 3);
		var r = 0.05 + 0.2 * random();
		c.y = r;
		var color = randomVec3Positive();
		scene.push(new Sphere(c, r, color));
	}


	var console = {
		timers: {},
		time: function(n) {
			this.timers[n] = performance.now();
		},
		timeEnd: function(n) {
			var ms = performance.now() - this.timers[n];
			this.log(n + " " + ms + " ms");
		},
		log: function() {
			var d = document.createElement('div');
			d.textContent = [].join.call(arguments, " ");
			window.debug.appendChild(d);
		}
	};


	function render() {
		var t = Date.now() / 3000.0
		camera.position.set(Math.cos(t), 0.5, Math.sin(t)).normalize().multiplyScalar(4.5);
		camera.target = new THREE.Vector3(0, 0.75, 0);
		camera.lookAt(camera.target);
		camera.updateProjectionMatrix();
		camera.updateMatrixWorld();

		window.debug.innerHTML = "";


		(console || window.console).time("voxelGrid build");

		var voxelGrid = new VoxelGrid(8, vec3(-1.1,-0.1,-1.1), vec3(2.2), 1, 8);
		bunnyTris.forEach(function(o) {
			voxelGrid.add(o);
		});

		(console || window.console).timeEnd("voxelGrid build");

		(console || window.console).time("voxelGrid2 build");

		var voxelGrid2 = new VoxelGrid(16, vec3(-8.1,-0.1,-8.1), vec3(16.2), 0, 1);
		scene.forEach(function(o) {
			voxelGrid2.add(o);
		});

		(console || window.console).timeEnd("voxelGrid2 build");


		console.time("trace");

		var rays = [];
		var epsilon = 0.0001;

		for (var y=0; y<canvasSize; y++) {
			for (var x=0; x<canvasSize; x++) {
				for (var dy=0; dy<AA_SIZE; dy++) {
					for (var dx=0; dx<AA_SIZE; dx++) {
						var r = setupRay(camera, x+dx/AA_SIZE, (canvasSize-y-1)+dy/AA_SIZE, canvasSize, canvasSize);
						rays.push(r);
					}
				}
			}
		}

		var plane = new Plane(vec3(0,0,0), vec3(0,1,0), vec3(0.5));
		var rayCount = rays.length;
		var lastRayCount = rayCount;
		VoxelGrid.stepCount = 0;
		VoxelGrid.cmpCount = 0;
		console.log("Tracing " + rays.length + " primary rays");
		for (var j=0; j<5; j++) {
			for (var i=0; i<rays.length; i++) {
				var r = rays[i];
				if (r.finished) continue;
				var hit = voxelGrid.intersect(r);	
				var hit1 = voxelGrid2.intersect(r);
				if (!hit || (hit1 && hit1.distance < hit.distance)) {
					hit = hit1;
				}
				var hit2 = plane.intersect(r);
				if (hit2 > 0 && (!hit || hit2 < hit.distance)) {
					hit = {distance: hit2, obj: plane};
				}
				if (hit) {
					if (r.bounce < 5) {
						rayCount++;
						r.o = add(r.o, mulS(r.d, hit.distance));
						var c = hit.obj.color(r);
						r.transmit = mul(r.transmit, c);
						var nml = hit.obj.normal(r.o);
						r.d = normalize(reflect(r.d, nml));
						r.o = add(r.o, mulS(r.d, epsilon));
						r.bounce++;
					} else {
						r.finished = true;
					}
				} else {
					var bg = mulS(vec3(0.4+sat(-r.d.y), 0.6, 0.8+sat(r.d.y)*abs(r.d.z)), 2);
					bg = add(bg, mulS(vec3(10.0, 6.0, 4.0), pow(sat(dot(r.d, normalize(vec3(5.0, 15.0, 10.0)))), 64.0) ));
					r.light = add(r.light, mul(r.transmit, bg));
					r.finished = true;
				}
			}
			console.log("Bounce " + j + ", traced " + (rayCount-lastRayCount));
			lastRayCount = rayCount;
		}

		console.timeEnd("trace");

		console.log("Traced " + rayCount + " rays");
		console.log(VoxelGrid.stepCount, "VoxelGrid steps");
		console.log(VoxelGrid.cmpCount, "Primitive intersection tests");
		console.log(VoxelGrid.stepCount / rays.length, "steps per ray")
		console.log(VoxelGrid.cmpCount / rays.length, "primitive intersections tests per ray");

		for (var i=0; i<canvasSize*canvasSize; i++) {
			var c = vec3();
			for (var dy=0; dy<AA_SIZE; dy++) {
				for (var dx=0; dx<AA_SIZE; dx++) {
					var r = rays[i*AA_SIZE*AA_SIZE+dy*AA_SIZE+dx];
					c = add(c, addS(neg(expV(neg(r.light))), 1));
				}
			}
			c = mulS(c, 1/(AA_SIZE*AA_SIZE));
			id.data[i*4 + 0] = c.x * 255;
			id.data[i*4 + 1] = c.y * 255;
			id.data[i*4 + 2] = c.z * 255;
			id.data[i*4 + 3] = 255;
		}

		ctx.putImageData(id, 0, 0);
	};

	var tick = function() {
		render();
		requestAnimationFrame(tick);
	};

	requestAnimationFrame(tick);
});

