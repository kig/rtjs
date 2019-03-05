
// Ray tracing setup

uniform float deviceEpsilonTrace;

uniform vec3 cameraFocusPoint;
uniform float focusDistance;
uniform float cameraApertureSize;

uniform float roughness;

uniform bool showFocalPlane;
uniform bool stripes;


const float SKY_DISTANCE = 1e6;


Ray setupRay(vec2 fragCoord) {
	vec2 uv = (fragCoord.xy / iResolution.xy) * 2.0 - 1.0;
	uv.x *= iResolution.x / iResolution.y;

	vec3 origin = cameraPosition;
	vec3 direction = normalize( unproject(vec3(uv, 1.0)) - origin );

	// Camera aperture simulation

	vec3 cameraFocusVector = cameraFocusPoint - origin;
	float focusDistance = dot(cameraFocusVector, direction);
	vec3 target = origin + (direction * focusDistance);

	origin += applyMatrix4(diskPoint(5.0+origin+direction) * cameraApertureSize, cameraMatrixWorld);
	direction = normalize(target - origin);

	return Ray(
		origin,
		direction,
		1.0 / direction,
		vec3(1.0),
		vec3(0.0, 0.0, 0.0),
		0.0,
		-1
	);
}

vec3 getColor(in Ray r, in int index) {
	if (index < 0) {
		return vec3(0.1);
	} else {
		return vec3(0.85, 0.53, 0.15);
	}
}

vec3 trace(Array vgArray, vec2 fragCoord) {
	float epsilon = deviceEpsilonTrace;
	Plane plane = Plane(vec3(0.0), vec3(0.0, 1.0, 0.0), vec3(0.5));
	Ray r = setupRay(fragCoord);

    int headOff = 18 * readInt(vgArray, 0) + 4;

	vec3 bg0;

	if (!costVis) {
		bg0 = vec3(0.6+clamp(-1.0, 0.0, 1.0), 0.7, 0.8+(0.4*0.0) * abs(0.0));
		bg0 += 0.25 * vec3(10.0, 6.0, 4.0) * 4.0 * pow(clamp(dot(vec3(0.0, 1.0, 0.0), normalize(vec3(6.0, 10.0, 8.0))), 0.0, 1.0), 64.0);
		bg0 += vec3(4.0, 5.0, 7.0) * abs(1.0 - 0.0);
		bg0 *= 0.2;
	}

	float bounce = 0.0;

	for (int j = 0; j < 6; j++) {
		Hit hit = setupHit();
		Hit hit2 = setupHit();
		// intersectSphere(r, vec3(0.0, 0.5, 0.0), 0.5, hit);
		intersectGridNode(vgArray, r, headOff, hit);
		intersectPlane(r, plane, hit2);
		if (hit2.distance < hit.distance) {
			hit = hit2;
		}
		if (hit.distance >= SKY_DISTANCE) {
			break;
		}
		if (showFocalPlane) {
			vec3 cameraFocusVector = cameraFocusPoint - r.o;
			float focusDistance = (dot(cameraFocusVector, r.d));
			r.light.r += 10.0 * (
				2.0 * clamp(0.1 / cameraApertureSize - abs(hit.distance - focusDistance), 0.0, 0.1) +
				200.0 * max(0.0, 0.03 - abs(hit.distance - focusDistance)) 
			);
		}
		bounce++;
		r.lastTested = hit.index;
		r.o = r.o + (r.d * hit.distance);
		r.pathLength += hit.distance;
		if (!costVis) {
			r.light += r.transmit * (1.0-exp(-hit.distance/40.0)) * bg0;
			vec3 c = getColor(r, hit.index);
			r.transmit = r.transmit * c;
		}
		vec3 nml = hit.index >= 0 ? triNormal(vgArray, r.o, hit.index) : plane.normal;
		// float troughness = mod(float(hit.index+2), 100.0) / 99.0; 
		// vec3 nml = hit.index >= 0 ? normalize(r.o-vec3(0.0, 0.5, 0.0)) : plane.normal;
		// r.d = normalize(reflect(r.d, nml));
		if (stripes) {
			r.d = normalize(mix(reflect(r.d, nml), nml + randomVec3(5.0+r.o+r.d), mod(roughness*dot(r.o, r.o)*10.0, 1.0)));// + (abs(dot(r.d, nml)) * roughness) * randomVec3(5.0+r.o+r.d));
		} else {
			r.d = normalize(mix(reflect(r.d, nml), nml + randomVec3(5.0+r.o+r.d), hit.index >= 0 ? roughness : 0.05));// + (abs(dot(r.d, nml)) * roughness) * randomVec3(5.0+r.o+r.d));
		}
		r.invD = 1.0 / r.d;
		r.o = r.o + nml * epsilon;
	}
if (!costVis) {
	if (bounce < 4.5) {
		vec3 bg = vec3(0.6+clamp(-r.d.y, 0.0, 1.0), 0.7, 0.8+(0.4*r.d.x) * abs(r.d.z));
		bg += vec3(10.0, 6.0, 4.0) * pow(clamp(dot(r.d, normalize(vec3(6.0, 10.0, 8.0))), 0.0, 1.0), 16.0);
		bg += vec3(10.0, 8.0, 6.0) * 4.0 * pow(clamp(dot(r.d, normalize(vec3(6.0, 10.0, 8.0))), 0.0, 1.0), 256.0);
		bg += vec3(4.0, 5.0, 7.0) * (abs(1.0 - r.d.z));
		r.light = r.light + r.transmit * bg;
	// } else {
		// r.light *= 10.0 * clamp((1.0-0.5*(length(r.o-cameraPosition)-3.5)), 0.0, 1.0);
	}
}
	return r.light;
}
