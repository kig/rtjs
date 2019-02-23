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

const float SKY_DISTANCE = 1e6;

struct Plane {
	vec3 point;
	vec3 normal;
	vec3 color;
};

float InterleavedGradientNoise(vec2 xy) {
  return fract(52.9829189f
              * fract(xy.x * 0.06711056f
                   + xy.y * 0.00583715f));
}

float random(float seed) {
	return InterleavedGradientNoise(gl_FragCoord.xy + seed);
}

vec3 randomVec3(vec3 p) {
	return vec3(random(p.x), random(p.y), random(p.z));
}

vec3 diskPoint(float seed) {
	float a = random(seed);
	float r = random(a);
	float sr = sqrt(r);
	return vec3(cos(a) * sr, sin(a) * sr, 0.0);
}

Hit setupHit() {
	return Hit(-1.0, 1.0e9);
}

Ray setupRay() {
	vec2 uv = (gl_FragCoord.xy / iResolution.xy) * 2.0 - 1.0;
	uv.x *= iResolution.x / iResolution.y;

	vec3 origin = cameraPosition;
	vec3 direction = normalize( (projectionMatrix * viewMatrix * vec4(vec3(uv.xy, 0.5) - origin, 1.0)).xyz );

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
		-1.0
	);
}

vec3 getColor(Ray r, float index) {
	if (index < 0.0) {
		return vec3(0.5);
	} else {
		return vec3(0.85, 0.53, 0.15);
	}
}

void intersectPlane(Ray ray, Plane p, inout Hit hit) {
	float pd = dot(p.normal, ray.d);
	if (abs(pd) > 0.00001) {
		float dist = dot(p.normal, p.point - ray.o) / pd;
		if (dist > 0.0) {
			hit.distance = dist;
		}
	}
}

vec3 trace(Array vgArray) {
	float epsilon = 0.0001;
	Plane plane = Plane(vec3(0.0), vec3(0.0, 1.0, 0.0), vec3(0.5));
	Ray r = setupRay();

	for (float j = 0.0; j < 6.0; j++) {
		Hit hit = setupHit();
		Hit hit2 = setupHit();
		//intersectGridLeaf(vgArray, r, 0.0, hit);
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
			vec3 nml = hit.index >= -1.0 ? triNormal(vgArray, r.o, hit.index) : plane.normal;
			r.d = normalize(reflect(r.d, nml) + randomVec3(r.o) * 0.1);
			r.o = r.o + nml * epsilon;
		} else {
			vec3 bg = vec3(0.6+clamp(-r.d.y, 0.0, 1.0), 0.7, 0.8+(0.4*r.d.x) * abs(r.d.z));
			bg += vec3(10.0, 6.0, 4.0) * 4.0 * pow(clamp(dot(r.d, normalize(vec3(6.0, 10.0, 8.0))), 0.0, 1.0), 64.0);
			//bg += vec3(3.0, 5.0, 7.0) * abs(1.0 - r.d.z);
			r.light = mix(r.light + (r.transmit * bg), bg, 1.0 - exp(-r.pathLength/40.0));
			break;
		}
	}
	return r.light;
}
