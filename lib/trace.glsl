
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

	origin += applyMatrix4(diskPoint(fragCoord) * cameraApertureSize, cameraMatrixWorld);
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
		return vec3(0.01);
	} else {
		return mod(r.o.z*-r.o.y*r.o.x, 0.05) > 0.005 ? vec3(0.95) : vec3(0.05);
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

    int headOff = 18 * readInt(vgArray, 0) + 4;

	vec3 bg0;

	if (!costVis) {
		bg0 = vec3(0.6+clamp(-1.0, 0.0, 1.0), 0.7, 0.8+(0.4*0.0) * abs(0.0));
		bg0 += 0.25 * vec3(10.0, 6.0, 4.0) * 4.0 * pow(clamp(dot(vec3(0.0, 1.0, 0.0), normalize(vec3(6.0, 10.0, 8.0))), 0.0, 1.0), 64.0);
		bg0 += vec3(4.0, 5.0, 7.0) * abs(1.0 - 0.0);
		bg0 *= 0.2;
	}

	float fresnel0;
	Hit hit0;
	vec3 nml0;

	vec3 basis[3];

	// TODO: Control number of rays at each bounce depth
	//
	// Sample primary rays for motion blur, AA, lens effects
	// Sample secondary rays for diffuse hemisphere light and reflections
	// Sample longer light paths for reflections, refractions, bounce light, SSS
	// Sample reverse light paths for caustics, difficult bounce light
	// Post-process for glow
	// More rays for less converged paths
	//
	// Trace light rays from emissive surfaces for BDPT.
	// Connect primary hit to surrounding secondary paths. (Covariance in MC)
	// Owen transformed Sobol sequence for faster numerical integration convergence.
	//
	// Find high energy transport path segments == path segment contribution to ray light.
	//   - Remove path segments one by one, reconnect path with shadow ray & pdf weight.
	//   - seg.contribution = (path.light - (path - seg).light) / path.light
	//   - Highest contributing segments are things like a ray through a door gap
	//   - Removing negative segments yields a shorter more energetic path

	bool hitScene = traceBounce(vgArray, headOff, r, plane, bg0, hit0, fresnel0, nml0, true);
	vec4 light = vec4(0.0);
	Ray primaryHit = r;

	float phi = 0.5 + 0.5 * sqrt(5.0);
	float rphi = 1.0 / phi;

	if (hitScene) {
		orthoBasis(basis, reflect(r.d, nml0));

		float fakeBounce = 10.0;
		float error = 1.0;
		vec3 offV;
		float roff = randomBlue((r.o+r.d+fakeBounce)).x * rphi;
		for (int j = 0; j < 4; j++) {
			// if (light.a > 4.0 && error < 0.01) {
			// 	// light = vec4(0.0, 10.0, 10.0, light.a);
			// 	break;
			// }
			if ((j & 3) == 0) offV = randomVec3(r.o+r.d+fakeBounce++);
			else if ((j & 1) == 0) offV = cross(offV, basis[2]);
			else if ((j & 3) == 3) offV = cross(offV, basis[0]);
			else cross(offV, basis[1]);
			offV = randomVec3(r.o+r.d+fakeBounce++);
			float fresnel = fresnel0;
			Hit hit = hit0;
			vec3 nml = nml0;
			for (int i = 0; i < 5; i++) {
				if (stripes) {
					r.d = normalize(mix(reflect(r.d, nml), nml + (i != 0 ? randomVec3(r.o+r.d+fakeBounce++) : offV), (1.0-fresnel) * mod(roughness*dot(r.o, r.o)*10.0, 1.0)));// + (abs(dot(r.d, nml)) * roughness) * randomVec3(5.0+r.o+r.d));
				} else {
					r.d = normalize(mix(reflect(r.d, nml), nml + (i != 0 ? randomVec3(r.o+r.d+fakeBounce++) : offV), (1.0-fresnel) * (hit.index >= 0 ? roughness : 0.05)));// + (abs(dot(r.d, nml)) * roughness) * randomVec3(5.0+r.o+r.d));
				}
				r.invD = 1.0 / r.d;
				r.o = r.o + nml * epsilon;
				if (!traceBounce(vgArray, headOff, r, plane, bg0, hit, fresnel, nml, false)) {
					vec3 bg = vec3(0.6+clamp(-r.d.y, 0.0, 1.0), 0.7, 0.8+(0.4*r.d.x) * abs(r.d.z));
					bg += vec3(10.0, 6.0, 4.0) * pow(clamp(dot(r.d, normalize(vec3(6.0, 10.0, 8.0))), 0.0, 1.0), 16.0);
					bg += vec3(10.0, 8.0, 6.0) * 4.0 * pow(clamp(dot(r.d, normalize(vec3(6.0, 10.0, 8.0))), 0.0, 1.0), 256.0);
					bg += vec3(4.0, 5.0, 7.0) * (abs(1.0 - r.d.z));
					vec3 rl = r.light + r.transmit * bg;
					vec3 ev = (light.rgb / max(1.0, light.a)) - ((light.rgb + rl) / (light.a + 1.0));
					error = dot(ev, ev);
					light.rgb += rl;
					light.a += 1.0;
					break;
				}
			}
			r = primaryHit;
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
