// Global config

uniform bool costVis;

uniform mat4 cameraInverseMatrix;
uniform mat4 cameraMatrixWorld;

uniform vec2 iResolution;
uniform float iTime;
uniform float iFrame;

// Phi (1 + sqrt 5) / 2
#define IR_1 1.618033988749895
// 1 + sqrt 2
#define IR_2 2.414213562373095
// (9 + sqrt 221) / 10
#define IR_3 2.386606874731851

#define PHI 1.618033988749895
#define RPHI 0.618033988749895

#define PI 3.141592653589793
#define E 2.718281828459045

#define SQRT1_2 0.7071067811865476
#define SQRT2 1.4142135623730951
#define SQRT3 1.7320508075688772
#define SQRT5 2.23606797749979

#define LN2 0.6931471805599453
#define LN10 2.302585092994046
#define LOG2E 1.4426950408889634
#define LOG10E 0.4342944819032518


// Structs

struct Hit {
    int index;
    float distance;
};

struct Ray {
    vec3 o;
    vec3 d;
    vec3 invD;
    vec3 transmit;
    vec3 light;
    float pathLength;
    float IOR;
    int lastTested;
};

struct Coating {
    // Specular material transmit
    vec3 color;

    // Specular material surface roughness
    float roughness;

    // Index of refraction of the specular material
    // P(specular) = pow((currentIOR-materialIOR)/(currentIOR+materialIOR), 2)
    float IOR;

    // Retroreflectiveness of the specular surface.
    // Retroreflect is a property of concave surface structures where 
    // the incoming light is reflected back to its source
    // after three+ bounces.
    float retroReflect;

    // Anisotropy in UV space
    // Negative values are used for circular pattern radius
    float anisotropy;
    float anisotropyAngle;

    // Density of the pigment for subsurface scattering
    float density;

    // Proportion of forward scattering
    float forwardScatter;

    // Thickness of the layer
    float thickness;

    // Proportion of light refracted through the material from surface
    float transmission;

    // Subsurface scattering probability per unit of depth
    float subsurface;
};

struct Material {
    Coating specular;
    Coating coat;

    Coating pigment;

    // Thin-shell object.
    float thinWalled;
    float wallThickness;

    // Color of light emitted by the material
    vec3 emission;
};

struct Array {
    int width;
};

struct Plane {
	vec3 point;
	vec3 normal;
	vec3 color;
};

Hit setupHit() {
	return Hit(-1, 1.0e9);
}


// Data access

uint readUint(usampler2D tex, int texWidth, int index) {
    int v = index / texWidth;
    int u = index - (v * texWidth);
    return texelFetch(tex, ivec2(u, v), 0).r;
}

int readInt(isampler2D tex, int texWidth, int index) {
    int v = index / texWidth;
    int u = index - (v * texWidth);
    return texelFetch(tex, ivec2(u, v), 0).r;
}

float readFloat(sampler2D tex, int texWidth, int index) {
    int v = index / texWidth;
    int u = index - (v * texWidth);
    return texelFetch(tex, ivec2(u, v), 0).r;
}

vec2 readVec2(sampler2D tex, int texWidth, int index) {
    return vec2(
        readFloat(tex, texWidth, index),
        readFloat(tex, texWidth, index + 1)
    );
}

vec3 readVec3(sampler2D tex, int texWidth, int index) {
    return vec3(
        readFloat(tex, texWidth, index),
        readFloat(tex, texWidth, index + 1),
        readFloat(tex, texWidth, index + 2)
    );
}


// Random numbers

uniform sampler2D blueNoise;

vec4 randomBlue(ivec2 st) {
	return texelFetch(blueNoise, st, 0);
}

float random(vec2 st) {
    return fract(fract(dot(sqrt(5.0) * (gl_FragCoord.xy / iResolution + st.xy + mod(vec2(iFrame/sqrt(2.0), iTime), vec2(1.0))), vec2(12.19898,78.233)))* 43758.5453123);
}

float InterleavedGradientNoise(in vec2 xy) {
  return fract(52.9829189f
              * fract(xy.x * 0.06711056f
                   + xy.y * 0.00583715f));
}

vec3 randomVec3U(in vec3 p) {
	float a = random(p.xz) * (2.0 * 3.14159);
	float r = random(p.yx);
    r = 2.0 * r - 1.0;
    return vec3(sqrt(1.0 - r * r) * vec2(cos(a), sin(a)), r);
}

vec3 diskPointU(in vec3 p) {
	float a = random(p.xz) * (2.0 * 3.14159);
	float r = random(p.yx);
	float sr = sqrt(r);
	return vec3(cos(a) * sr, sin(a) * sr, 0.0);
}

vec3 randomVec3(in ivec2 st) {
	vec2 ar = randomBlue(st).zw;
	float a = ar.x * (2.0 * 3.14159);
	float r = ar.y;
    r = 2.0 * r - 1.0;
    return vec3(sqrt(1.0 - r * r) * vec2(cos(a), sin(a)), r);
}

vec3 randomInUnitSphere(in ivec2 st) {
	vec3 ar = randomBlue(st).zwx;
	float a = ar.x * (2.0 * 3.14159);
	float r = ar.y;
    r = 2.0 * r - 1.0;
    return vec3(sqrt(1.0 - r * r) * vec2(cos(a), sin(a)), r) * pow(ar.z, 1.0 / 3.0);
}

vec3 randomVec3(in vec3 p) {
	return randomVec3(
        ivec2(
            mod(sqrt(5.0) * (
                gl_FragCoord.xy + 
                43758.5453123 * dot(vec2(12.9898,78.233), p.xz) + 
                vec2(1819278.233 * p.y + iFrame, iTime)),

                vec2(1024.0)
            )
        )
    );
}

vec4 randomBlue(in vec3 p) {
    return randomBlue(ivec2(mod(gl_FragCoord.xy + 43758.5453123*vec2(12.9898,78.233)*p.xz + vec2(1819278.233*p.y+iFrame, floor(iFrame/1024.0)), vec2(1024.0))));
}

vec3 diskPoint(in vec2 p) {
	vec2 ar = randomBlue(ivec2(mod(43758.5453123 * p * 4.0 + vec2(12.9898,78.233)*vec2(iFrame, floor(iFrame/1024.0)), vec2(1024.0)))).xy;
	ar.x *= (2.0 * 3.14159);
	float sr = sqrt(ar.y);
	return vec3(cos(ar.x) * sr, sin(ar.x) * sr, 0.0);
}

// vec3 randomVec3(in vec3 p) {
//     if (iFrame > 30.0) {
//         return randomVec3U(p);
//     } else {
//         return randomVec3B(p);
//     }
// }

// Three.js math

vec3 applyMatrix4(in vec3 v, in mat4 m) {
	float w = 1.0 / ( m[0][3] * v.x + m[1][3] * v.y + m[2][3] * v.z + m[3][3] );

	return vec3(
		( m[0][0] * v.x + m[1][0] * v.y + m[2][0] * v.z + m[3][0] ) * w,
		( m[0][1] * v.x + m[1][1] * v.y + m[2][1] * v.z + m[3][1] ) * w,
		( m[0][2] * v.x + m[1][2] * v.y + m[2][2] * v.z + m[3][2] ) * w
	);
}

vec3 applyMatrix4Rot(in vec3 v, in mat4 m) {
	float w = 1.0; // / ( m[0][3] * v.x + m[1][3] * v.y + m[2][3] * v.z + m[3][3] );

	return vec3(
		( m[0][0] * v.x + m[1][0] * v.y + m[2][0] * v.z ) * w,
		( m[0][1] * v.x + m[1][1] * v.y + m[2][1] * v.z ) * w,
		( m[0][2] * v.x + m[1][2] * v.y + m[2][2] * v.z ) * w
	);
}

vec3 unproject(in vec3 v) {
	return applyMatrix4(v, cameraInverseMatrix);
}

void orthoBasis(out vec3 basis[3], vec3 n)
{
	basis[2] = vec3(n.x, n.y, n.z);
	basis[1] = vec3(0.0, 0.0, 0.0);

	if ((n.x < 0.6) && (n.x > -0.6))
		basis[1].x = 1.0;
	else if ((n.y < 0.6) && (n.y > -0.6))
		basis[1].y = 1.0;
	else if ((n.z < 0.6) && (n.z > -0.6))
		basis[1].z = 1.0;
	else
		basis[1].x = 1.0;


	basis[0] = cross(basis[1], basis[2]);
	basis[0] = normalize(basis[0]);

	basis[1] = cross(basis[2], basis[0]);
	basis[1] = normalize(basis[1]);

}

vec2 toOct(in vec3 v) {
	vec3 p = v * (1.0 / (abs(v.x) + abs(v.y) + abs(v.z)));
	if (p.z < 0.0) {
		float px = p.x;
		p.x = (1.0 - abs(p.y)) * (p.x >= 0.0 ? 1.0 : -1.0);
		p.y = (1.0 - abs(px)) * (p.y >= 0.0 ? 1.0 : -1.0);
	}
    return p.xy * 0.5 + 0.5;
}

vec3 fromOct(in vec2 uv) {
    uv = uv * 2.0 - 1.0;
	vec3 p = vec3(uv, 1.0 - abs(uv.x) - abs(uv.y));
	if (p.z < 0.0) {
		float px = p.x;
		p.x = (1.0 - abs(p.y)) * (p.x >= 0.0 ? 1.0 : -1.0);
		p.y = (1.0 - abs(px)) * (p.y >= 0.0 ? 1.0 : -1.0);
	}
	return normalize(p);
}

vec3 fromGrid(vec3 coord, float scale, vec3 origin) {
    return origin + (coord * scale);
}

vec3 toGrid(vec3 point, float scale, vec3 origin) {
    return (point - origin) / scale;
}


// Ray tracing primitives

uniform float deviceEpsilon;

void intersectPlane(in Ray ray, in Plane p, inout Hit hit) {
	float pd = dot(p.normal, ray.d);
	if (abs(pd) > 0.00001) {
		float dist = dot(p.normal, p.point - ray.o) / pd;
		if (dist > 0.0 && length(p.point - (ray.o + ray.d * dist)) < 2.0) {
			hit.distance = dist;
			hit.index = -2;
		}
	}
}

void intersectSphere(in Ray ray, in vec3 p, in float r, inout Hit hit) {
	vec3 rc = ray.o - p; 
	float c = dot(rc, rc) - r*r;
	float b = dot(ray.d, rc);
	float d = b*b - c;
	if (d < 0.0) {
		return;
	}
	float t = -b - sqrt(d);
	if (t <= 0.0) {
		return;
	}
	hit.distance = t;
	hit.index = -3;
}

float forceIntersectBox(in Ray ray, in vec3 origin, in vec3 dims) {
    vec3 bmin = origin;
    vec3 bmax = origin + dims;
    vec3 t0 = (bmin - ray.o) * ray.invD;
    vec3 t1 = (bmax - ray.o) * ray.invD;
    vec3 swaps = (sign(ray.invD) + 1.0) * 0.5;
    vec3 swaps2 = 1.0 - swaps;
    vec3 t0_ = t0 * swaps + t1 * swaps2;

    return max(0.0, max(t0_.x, max(t0_.y, t0_.z)));
}

float intersectBox(in Ray ray, in vec3 origin, in vec3 dims) {
    vec3 bmin = origin;
    vec3 bmax = origin + dims;
    vec3 t0 = (bmin - ray.o) * ray.invD;
    vec3 t1 = (bmax - ray.o) * ray.invD;
    vec3 swaps = (sign(ray.invD) + 1.0) * 0.5;
    vec3 swaps2 = 1.0 - swaps;
    vec3 t0_ = t0 * swaps + t1 * swaps2;
    vec3 t1_ = t0 * swaps2 + t1 * swaps;

    float tmin = max(t0_.x, max(t0_.y, t0_.z));
    float tmax = min(t1_.x, min(t1_.y, t1_.z));

    if (tmax <= tmin) {
        return -1.0;
    }
    return max(0.0, tmin);
}

void intersectTri(vec3 v0, vec3 e1, vec3 e2, in int triIndex, inout Ray ray, inout Hit closestHit) {
    vec3 h = cross(ray.d, e2);
    float a = dot(e1, h);

    if (a > -0.00001 && a < 0.00001) {
        return;
    }

    float f = 1.0 / a;
    vec3 s = ray.o - v0;
    float u = f * dot(s, h);

    if (u < 0.0 - deviceEpsilon || u > 1.0 + deviceEpsilon) {
        return;
    }

    vec3 q = cross(s, e1);
    float v = f * dot(ray.d, q);

    if (v < 0.0 - deviceEpsilon || u + v > 1.0 + deviceEpsilon*2.0) {
        return;
    }

    // at this stage we can compute t to find out where
    // the intersection point is on the line
    float t = f * dot(e2, q);

    if (t < 0.000001 || t > closestHit.distance - deviceEpsilon*0.1) {
        return;
    }

    // if (costVis) {
    //     float borderWidth = 0.025;
    //     if (u < borderWidth || v < borderWidth || u+v > 1.0-borderWidth) {
    //         ray.light.g += 1.0;
    //         ray.light.b += 1.0;
    //     }
    // }

    closestHit.index = triIndex;
    closestHit.distance = t;
}

void intersectTri(in sampler2D triangles, in int trianglesWidth, inout Ray ray, in int triIndex, inout Hit closestHit) {
    int off = triIndex * 9;
    vec3 v0 = readVec3(triangles, trianglesWidth, off);
    vec3 e1 = readVec3(triangles, trianglesWidth, off+3) - v0;
    vec3 e2 = readVec3(triangles, trianglesWidth, off+6) - v0;

    intersectTri(v0, e1, e2, triIndex, ray, closestHit);
}

vec3 triNormal(
    in sampler2D triangles, in int trianglesWidth, 
    in sampler2D normals, in int normalsWidth,
    in vec3 point, in int triIndex
) {
    int off = triIndex * 9;
    int nmlOff = triIndex * 6;

    vec3 v0 = readVec3(triangles, trianglesWidth, off);
    vec3 e1 = readVec3(triangles, trianglesWidth, off+3) - v0;
    vec3 e2 = readVec3(triangles, trianglesWidth, off+6) - v0;
    
    float u = clamp(dot(point, e1), 0.0, 1.0);
    float v = clamp(dot(point, e2), 0.0, 1.0);

    vec3 n0 = fromOct(readVec2(normals, normalsWidth, nmlOff));
    vec3 n1 = fromOct(readVec2(normals, normalsWidth, nmlOff + 2));
    vec3 n2 = fromOct(readVec2(normals, normalsWidth, nmlOff + 4));

    return normalize(mix(mix(n0, n2, v), n1, u));
}

uint readMaterialIndex(in usampler2D tex, in int texWidth, in int idx) {
    return readUint(tex, texWidth, idx);
}

Coating newCoating() {
    return Coating(
        vec3(0.),
        0.,0.,0.,0.,
        0.,0.,0.,0.,
        0., 0.
    );
}

Material newMaterial() {
    Coating c = newCoating();
    return Material(c,c,c, 0.,0., vec3(0.));
}

Coating readCoating(in sampler2D tex, in int materialsWidth, inout int idx) {
    Coating c = newCoating();
    c.color = readVec3(tex, materialsWidth, idx); // 3
    idx += 3;
    c.roughness = readFloat(tex, materialsWidth, idx++); // 4
    c.IOR = readFloat(tex, materialsWidth, idx++);
    c.retroReflect = readFloat(tex, materialsWidth, idx++); // 6
    c.anisotropy = readFloat(tex, materialsWidth, idx++);
    c.anisotropyAngle = readFloat(tex, materialsWidth, idx++); // 8
    c.density = readFloat(tex, materialsWidth, idx++);
    c.forwardScatter = readFloat(tex, materialsWidth, idx++); // 10
    c.thickness = readFloat(tex, materialsWidth, idx++);
    c.transmission = readFloat(tex, materialsWidth, idx++); // 12
    c.subsurface = readFloat(tex, materialsWidth, idx++); // 13
    return c;
}

Material readMaterial(in sampler2D tex, in int materialsWidth, in int idx) {
    idx = idx * 44;
    Material m = newMaterial();
    m.specular = readCoating(tex, materialsWidth, idx); // 13
    m.coat = readCoating(tex, materialsWidth, idx); // 26
    m.pigment = readCoating(tex, materialsWidth, idx); // 39
    m.thinWalled = readFloat(tex, materialsWidth, idx++); // 40
    m.wallThickness = readFloat(tex, materialsWidth, idx++); // 41
    m.emission = readVec3(tex, materialsWidth, idx); // 44
    idx += 3;
    return m;
}