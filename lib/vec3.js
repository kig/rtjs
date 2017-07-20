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
	return v < 0 ? 0 : (v > 1 ? 1 : v);
};

var satV = function(v) {
	return vec3(sat(v.x), sat(v.y), sat(v.z));
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