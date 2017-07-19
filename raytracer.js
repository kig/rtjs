var vec3 = function(x, y, z) {
	x = x === undefined ? 0 : x;
	y = y === undefined ? x : y;
	z = z === undefined ? y : z;
	return {x: x, y: y, z: z};
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

var add = function(u,v) {
	return vec3(u.x+v.x, u.y+v.y, u.z+v.z);
};

var sub = function(u,v) {
	return vec3(u.x-v.x, u.y-v.y, u.z-v.z);
};

var mulS = function(u,s) {
	return vec3(u.x*s, u.y*s, u.z*s);
};

var addS = function(u,s) {
	return vec3(u.x+s, u.y+s, u.z+s);
};

var subS = function(u,s) {
	return vec3(u.x-s, u.y-s, u.z-s);
};

var dot = function(u,v) {
	return (u.x*v.x + u.y*v.y + u.z*v.z);
};

var reflect = function(v, nml) {
	return add(v, mulS(nml, 2));
};

var length = function(v) {
	return sqrt(dot(v,v));
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

var floorV = function(v) {
	return vec3(floor(v.x), floor(v.y), floor(v.z));
};

var ceilV = function(v) {
	return vec3(ceil(v.x), ceil(v.y), ceil(v.z));
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
	return {
		min: sub(this.center, vec3(this.radius, this.radius, this.radius)),
		max: add(this.center, vec3(this.radius, this.radius, this.radius)),
	};
};

Sphere.prototype.intersectBox = function(box) {
	return true;
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
	return {
		min: vec3(
			min(this._vertices[0].x, this._vertices[1].x, this._vertices[2].x),
			min(this._vertices[0].y, this._vertices[1].y, this._vertices[2].y),
			min(this._vertices[0].z, this._vertices[1].z, this._vertices[2].z)
			),
		max: vec3(
			max(this._vertices[0].x, this._vertices[1].x, this._vertices[2].x),
			max(this._vertices[0].y, this._vertices[1].y, this._vertices[2].y),
			max(this._vertices[0].z, this._vertices[1].z, this._vertices[2].z)
			)
	};
};

Triangle.prototype.intersectBox = function(box) {
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
						v = this._voxels[z][y][x] = this._leaf === 0 ? new Voxel() : new VoxelGrid(16, this.fromGrid(vec3(x,y,z)), vec3(this._scale), this._leaf-1, this._levelSize);
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

	for (var i = 0; i < verts.length; i += 3*3) {
		var u = vec3(verts[i], verts[i+1], verts[i+2]);
		var v = vec3(verts[i+3], verts[i+4], verts[i+5]);
		var w = vec3(verts[i+6], verts[i+7], verts[i+8]);
		var x = vec3(normals[i], normals[i+1], normals[i+2]);
		var y = vec3(normals[i+3], normals[i+4], normals[i+5]);
		var z = vec3(normals[i+6], normals[i+7], normals[i+8]);
		scene.push(new Triangle([u,v,w], [x,y,z], add(vec3(0.8), mulS(randomVec3Positive(), 0.2))));
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
		camera.position.set(Math.cos(t), 0.5, Math.sin(t)).normalize().multiplyScalar(4.0);
		camera.target = new THREE.Vector3(0, 0.5, 0);
		camera.lookAt(camera.target);
		camera.updateProjectionMatrix();
		camera.updateMatrixWorld();

		window.debug.innerHTML = "";

		console.time("voxelGrid build");

		var voxelGrid = new VoxelGrid(11, vec3(-8), vec3(16), 1, 11);
		scene.forEach(function(o) {
			voxelGrid.add(o);
		});

		console.timeEnd("voxelGrid build");


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
		var lastRayCount = 0;
		VoxelGrid.stepCount = 0;
		VoxelGrid.cmpCount = 0;
		console.log("Tracing " + rays.length + " primary rays");
		for (var j=0; j<5; j++) {
			for (var i=0; i<rays.length; i++) {
				var r = rays[i];
				if (r.finished) continue;
				var hit = voxelGrid.intersect(r);
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
					c.x += (1.0-exp(-r.light.x));
					c.y += (1.0-exp(-r.light.y));
					c.z += (1.0-exp(-r.light.z));
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

