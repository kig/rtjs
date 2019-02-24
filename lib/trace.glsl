// uniform mat4 modelMatrix; // = object.matrixWorld
// uniform mat4 modelViewMatrix; // = camera.matrixWorldInverse * object.matrixWorld
uniform mat4 projectionMatrix; // = camera.projectionMatrix
// uniform mat4 viewMatrix; // = camera.matrixWorldInverse
// uniform mat3 normalMatrix; // = inverse transpose of modelViewMatrix
// uniform vec3 cameraPosition; // = camera position in world space

uniform vec3 cameraFocusPoint;
uniform float focusDistance;
uniform float cameraApertureSize;

uniform vec2 iResolution;

uniform float iTime;

const float SKY_DISTANCE = 1e6;

struct Plane {
	vec3 point;
	vec3 normal;
	vec3 color;
};

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
	float a = random(p.xz);
	float r = random(p.yx);
	float h = random(vec2(a,r)) * 2.0 - 1.0;
	float sr = sqrt(r);
	return vec3(cos(a) * sr, sin(a) * sr, h);
}

// vec3 diskPoint(in float seed) {
// 	float a = random(seed);
// 	float r = random(a);
// 	float sr = sqrt(r);
// 	return vec3(cos(a) * sr, sin(a) * sr, 0.0);
// }

Hit setupHit() {
	return Hit(-1, 1.0e9);
}

vec3 unproject(in vec3 v) {
	vec4 v1 = vec4(v, 1.0) * projectionMatrix;
	vec4 v2 = v1 * viewMatrix;
	return v2.xyz;
}

Ray setupRay(vec2 fragCoord) {
	vec2 uv = (fragCoord.xy / iResolution.xy) * 2.0 - 1.0;
	uv.x *= iResolution.x / iResolution.y;

	vec3 origin = cameraPosition;
	// vec3 direction = normalize (vec3(uv, 1.0) - cameraPosition);
	vec3 direction = normalize( unproject(vec3(uv, 0.5)) - origin );

	// Camera aperture simulation

	// float focusDistance = length(origin - cameraFocusPoint);
	// vec3 target = origin + (direction * focusDistance);

	// origin += (diskPoint() * cameraApertureSize) * cameraMatrixWorld;
	// direction = normalize(target - origin); 

	return Ray(
		origin,
		direction,
		vec3(1.0),
		vec3(0.0),
		0.0,
		-1
	);
}

vec3 getColor(in Ray r, in int index) {
	if (index < 0) {
		return vec3(0.5);
	} else {
		return vec3(0.85, 0.53, 0.15);
	}
}

void intersectPlane(in Ray ray, in Plane p, inout Hit hit) {
	float pd = dot(p.normal, ray.d);
	if (abs(pd) > 0.00001) {
		float dist = dot(p.normal, p.point - ray.o) / pd;
		if (dist > 0.0) {
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

vec3 trace(Array vgArray, vec2 fragCoord) {
	float epsilon = 0.0001;
	Plane plane = Plane(vec3(0.0), vec3(0.0, 1.0, 0.0), vec3(0.5));
	Ray r = setupRay(fragCoord);

    int headOff = 18 * readInt(vgArray, 0) + 4;

	for (float j = 0.0; j < 6.0; j++) {
		Hit hit = setupHit();
		Hit hit2 = setupHit();
		// intersectSphere(r, vec3(0.0, 0.5, 0.0), 0.5, hit);
		intersectGridLeaf(vgArray, r, headOff, hit);
		intersectPlane(r, plane, hit2);
		if (hit2.distance < hit.distance) {
			hit = hit2;
		}
		if (hit.distance < SKY_DISTANCE) {
			r.lastTested = hit.index;
			r.o = r.o + (r.d * hit.distance);
			r.pathLength += hit.distance;
			vec3 c = getColor(r, hit.index);
			r.transmit = r.transmit * c;
			vec3 nml = hit.index >= 0 ? triNormal(vgArray, r.o, hit.index) : plane.normal;
			// vec3 nml = hit.index >= 0 ? normalize(r.o-vec3(0.0, 0.5, 0.0)) : plane.normal;
			r.d = normalize(reflect(r.d, nml) + (abs(dot(r.d, nml)) * 0.1) * randomVec3(r.o+r.d));
			r.o = r.o + nml * epsilon;
		} else {
			vec3 bg = vec3(0.6+clamp(-r.d.y, 0.0, 1.0), 0.7, 0.8+(0.4*r.d.x) * abs(r.d.z));
			bg += vec3(10.0, 6.0, 4.0) * 4.0 * pow(clamp(dot(r.d, normalize(vec3(6.0, 10.0, 8.0))), 0.0, 1.0), 64.0);
			bg += vec3(3.0, 5.0, 7.0) * abs(1.0 - r.d.z);
			r.light = mix(r.light + (r.transmit * bg), bg, 1.0 - exp(-r.pathLength/40.0));
			break;
		}
	}
	return r.light;
}
