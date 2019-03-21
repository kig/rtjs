
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

	vec3 dp = diskPoint(fragCoord) * cameraApertureSize;
	float dx = dp.x * dp.x;
	float dy = dp.y * dp.y;

	vec3 side = unproject(vec3(1.0, 0.0, 0.0));
	vec3 up = unproject(vec3(0.0, 1.0, 0.0));

	origin += applyMatrix4(dp, cameraMatrixWorld);
	// direction = normalize(target-origin); 
	direction = normalize(mix(mix(target - origin, sign(dp.x) * side, pow(dx, 8.0)), sign(dp.y) * up, pow(dy, 8.0)));

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
		return vec3(0.01);
	} else {
		return mix(vec3(0.95, 0.66, 0.15), vec3(0.05), float(stripes) * pow(fract(roughness*dot(r.o, r.o)*10.0), 0.125));
	}
}

bool traceBounce(in Array vgArray, in int headOff, inout Ray r, in Plane plane, in vec3 bg0, out Hit hit, out float fresnel, out vec3 nml, in bool isPrimaryRay) {
	hit = setupHit();
	Hit hit2 = setupHit();
	// intersectSphere(r, vec3(0.0, 0.5, 0.0), 0.5, hit);
	// intersectGrid(r, hit);
	intersectGridNode(vgArray, r, headOff, hit);
	intersectPlane(r, plane, hit2);
	if (hit2.distance < hit.distance) {
		hit = hit2;
	}
	if (hit.distance >= SKY_DISTANCE) {
		return false;
	}
	if (showFocalPlane && isPrimaryRay) {
		vec3 cameraFocusVector = cameraFocusPoint - r.o;
		float focusDistance = (dot(cameraFocusVector, r.d));
		r.light.r += 10.0 * (
			2.0 * clamp(0.1 / cameraApertureSize - abs(hit.distance - focusDistance), 0.0, 0.1) +
			200.0 * max(0.0, 0.03 - abs(hit.distance - focusDistance)) 
		);
	}
	r.lastTested = hit.index;
	r.o = r.o + (r.d * hit.distance);
	r.pathLength += hit.distance;
	nml = hit.index >= 0 ? triNormal(vgArray, r.o, hit.index) : plane.normal;
	// nml = hit.index >= 0 ? triNormal(triangles, trianglesWidth, normals, normalsWidth, r.o, hit.index) : plane.normal;
	// float troughness = mod(float(hit.index+2), 100.0) / 99.0; 
	// vec3 nml = hit.index >= 0 ? normalize(r.o-vec3(0.0, 0.5, 0.0)) : plane.normal;
	// r.d = normalize(reflect(r.d, nml));
	fresnel = pow(1.0 - abs(dot(r.d, nml)), 5.0);
	if (!costVis) {
		r.light += r.transmit * (1.0-exp(-hit.distance/40.0)) * bg0;
		vec3 c = getColor(r, hit.index);
		vec3 filmColor = abs(sin(r.o+4.0*r.d));
		r.transmit = r.transmit * c; //mix(c, filmColor, fresnel);
	}
	return true;
}

vec3 trace(Array vgArray, vec2 fragCoord) {
	float epsilon = deviceEpsilonTrace;

	Plane plane = Plane(vec3(0.0), vec3(0.0, 1.0, 0.0), vec3(0.5));
	Ray r = setupRay(fragCoord);

	vec4 light = vec4(0.0);

	int headOff = 18 * readInt(vgArray, 0) + 4;

	vec3 bg0;

	if (!costVis) {
		bg0 = vec3(0.6+clamp(-1.0, 0.0, 1.0), 0.7, 0.8+(0.4*0.0) * abs(0.0));
		bg0 += 0.25 * vec3(10.0, 6.0, 4.0) * 4.0 * pow(clamp(dot(vec3(0.0, 1.0, 0.0), normalize(vec3(6.0, 10.0, 8.0))), 0.0, 1.0), 64.0);
		bg0 += vec3(4.0, 5.0, 7.0) * abs(1.0 - 0.0);
		bg0 *= 0.2;
	}

	float fresnel;
	Hit hit;
	vec3 nml;

	bool hitScene = traceBounce(vgArray, headOff, r, plane, bg0, hit, fresnel, nml, true);

	if (hitScene) {
		float fakeBounce = 10.0;
		for (int i = 0; i < 5; i++) {
			if (stripes) {
				r.d = normalize(mix(reflect(r.d, nml), nml + randomVec3(r.o+r.d+fakeBounce++), (1.0-fresnel) * mod(roughness*dot(r.o, r.o)*10.0, 1.0)));// + (abs(dot(r.d, nml)) * roughness) * randomVec3(5.0+r.o+r.d));
			} else {
				r.d = normalize(mix(reflect(r.d, nml), nml + randomVec3(r.o+r.d+fakeBounce++), (1.0-fresnel) * (hit.index >= 0 ? roughness : 0.05)));// + (abs(dot(r.d, nml)) * roughness) * randomVec3(5.0+r.o+r.d));
			}
			r.invD = 1.0 / r.d;
			r.o = r.o + nml * epsilon;
			if (!traceBounce(vgArray, headOff, r, plane, bg0, hit, fresnel, nml, false)) {
				vec3 bg = vec3(0.6+clamp(-r.d.y, 0.0, 1.0), 0.7, 0.8+(0.4*r.d.x) * abs(r.d.z));
				bg += vec3(10.0, 6.0, 4.0) * pow(clamp(dot(r.d, normalize(vec3(6.0, 10.0, 8.0))), 0.0, 1.0), 16.0);
				bg += vec3(10.0, 8.0, 6.0) * 4.0 * pow(clamp(dot(r.d, normalize(vec3(6.0, 10.0, 8.0))), 0.0, 1.0), 256.0);
				bg += vec3(4.0, 5.0, 7.0) * (abs(1.0 - r.d.z));
				light.rgb += r.light + r.transmit * bg;
				light.a += 1.0;
				break;
			}
		}
	} else {
		vec3 bg = vec3(0.6+clamp(-r.d.y, 0.0, 1.0), 0.7, 0.8+(0.4*r.d.x) * abs(r.d.z));
		bg += vec3(10.0, 6.0, 4.0) * pow(clamp(dot(r.d, normalize(vec3(6.0, 10.0, 8.0))), 0.0, 1.0), 16.0);
		bg += vec3(10.0, 8.0, 6.0) * 4.0 * pow(clamp(dot(r.d, normalize(vec3(6.0, 10.0, 8.0))), 0.0, 1.0), 256.0);
		bg += vec3(4.0, 5.0, 7.0) * (abs(1.0 - r.d.z));
		light.rgb += r.light + r.transmit * bg;
		light.a += 1.0;
	}

	return light.rgb / max(1.0, light.a);
}
