// Global config

uniform bool costVis;

uniform mat4 cameraInverseMatrix;
uniform mat4 cameraMatrixWorld;

uniform vec2 iResolution;
uniform float iTime;

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
    int lastTested;
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

uniform sampler2D arrayTex;
uniform highp isampler2D iarrayTex;
uniform highp int arrayTexWidth;

int readInt(Array array, int index) {
    int v = index / array.width;
    int u = index - (v * array.width);
    return texelFetch(iarrayTex, ivec2(u, v), 0).r;
}

float readFloat(Array array, int index) {
    int v = index / array.width;
    int u = index  - (v * array.width);
    return texelFetch(arrayTex, ivec2(u, v), 0).r;
}

vec3 readVec3(Array array, int index) {
    return vec3(
        readFloat(array, index),
        readFloat(array, index + 1),
        readFloat(array, index + 2)
    );
}


// Random numbers

float InterleavedGradientNoise(in vec2 xy) {
  return fract(52.9829189f
              * fract(xy.x * 0.06711056f
                   + xy.y * 0.00583715f));
}

float random(vec2 st) {
    return fract(fract(dot(gl_FragCoord.xy/iResolution.xy + st.xy,
                         vec2(12.9898,78.233)))*
        43758.5453123);
}

vec3 randomVec3(in vec3 p) {
	float a = random(p.xz) * (2.0 * 3.14159);
	float r = random(p.yx);
    r = 2.0 * r - 1.0;
    return vec3(sqrt(1.0 - r * r) * vec2(cos(a), sin(a)), r);
}

vec3 diskPoint(in vec3 p) {
	float a = random(p.xz) * (2.0 * 3.14159);
	float r = random(p.yx);
	float sr = sqrt(r);
	return vec3(cos(a) * sr, sin(a) * sr, 0.0);
}


// Three.js math

vec3 applyMatrix4(in vec3 v, in mat4 m) {
	float w = 1.0 / ( m[0][3] * v.x + m[1][3] * v.y + m[2][3] * v.z + m[3][3] );

	return vec3(
		( m[0][0] * v.x + m[1][0] * v.y + m[2][0] * v.z + m[3][0] ) * w,
		( m[0][1] * v.x + m[1][1] * v.y + m[2][1] * v.z + m[3][1] ) * w,
		( m[0][2] * v.x + m[1][2] * v.y + m[2][2] * v.z + m[3][2] ) * w
	);
}

vec3 unproject(in vec3 v) {
	return applyMatrix4(v, cameraInverseMatrix);
}


// Ray tracing primitives

uniform float deviceEpsilon;

void intersectPlane(in Ray ray, in Plane p, inout Hit hit) {
	float pd = dot(p.normal, ray.d);
	if (abs(pd) > 0.00001) {
		float dist = dot(p.normal, p.point - ray.o) / pd;
		if (dist > 0.0 && length(ray.o + ray.d * dist) < 1.5) {
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

void intersectTri(in Array array, inout Ray ray, in int triIndex, inout Hit closestHit) {
        int off = 4 + triIndex * 18;
        vec3 v0 = readVec3(array, off);
        vec3 e1 = readVec3(array, off+3);
        vec3 e2 = readVec3(array, off+6);

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

        if (costVis) {
            float borderWidth = 0.025;
            if (u < borderWidth || v < borderWidth || u+v > 1.0-borderWidth) {
                ray.light.g += 1.0;
                ray.light.b += 1.0;
            }
        }

        closestHit.index = triIndex;
        closestHit.distance = t;
}

vec3 triNormal(in Array array, in vec3 point, in int triIndex) {
    int off = 4 + triIndex * 18;

    vec3 e1 = readVec3(array, off+3);
    vec3 e2 = readVec3(array, off+6);
    
    float u = clamp(dot(point, e1), 0.0, 1.0);
    float v = clamp(dot(point, e2), 0.0, 1.0);
    vec3 n0 = readVec3(array, off+9);
    vec3 n1 = readVec3(array, off+12);
    vec3 n2 = readVec3(array, off+15);

    return normalize(mix(mix(n0, n2, v), n1, u));
}
