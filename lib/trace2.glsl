
// Ray tracing setup

uniform float deviceEpsilonTrace;

uniform float cameraFocusDistance;
uniform float focusDistance;
uniform float cameraApertureSize;

uniform float roughness;

uniform bool showFocalPlane;
uniform bool stripes;

const float SKY_DISTANCE = 1e6;


Ray setupRay(vec2 fragCoord, float off) {
	vec2 uv = (fragCoord.xy / iResolution.xy) * 2.0 - 1.0;
	uv.x *= iResolution.x / iResolution.y;

	vec3 origin = cameraPosition;
	vec3 direction = normalize( unproject(vec3(uv, 1.0)) - origin );

	// Camera aperture simulation

	vec3 target = origin + (direction * cameraFocusDistance);

	origin += applyMatrix4(diskPoint(fragCoord+off) * cameraApertureSize, cameraMatrixWorld);
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
		return mix(vec3(0.95, 0.66, 0.15), vec3(0.05), float(stripes) * pow(fract(roughness*dot(r.o, r.o)*10.0), 0.125));
	}
}

bool traceBounce(inout Ray r, in Plane plane, in vec3 bg0, out Hit hit, out float fresnel, out vec3 nml, in bool isPrimaryRay) {
	hit = setupHit();
	Hit hit2 = setupHit();
	// intersectSphere(r, vec3(0.0, 0.5, 0.0), 0.5, hit);
	intersectGrid(r, hit);
	intersectPlane(r, plane, hit2);
	if (hit2.distance < hit.distance) {
		hit = hit2;
	}
	if (hit.distance >= SKY_DISTANCE) {
		return false;
	}
	if (showFocalPlane && isPrimaryRay) {
		r.light.r += 10.0 * (
			2.0 * clamp(0.1 / cameraApertureSize - abs(hit.distance - cameraFocusDistance), 0.0, 0.1) +
			200.0 * max(0.0, 0.03 - abs(hit.distance - cameraFocusDistance)) 
		);
	}
	r.lastTested = hit.index;
	r.o = r.o + (r.d * hit.distance);
	r.pathLength += hit.distance;
	nml = hit.index >= 0 ? triNormal(triangles, trianglesWidth, normals, normalsWidth, r.o, hit.index) : plane.normal;
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

vec3 trace(vec2 fragCoord) {
	float epsilon = deviceEpsilonTrace;

	Plane plane = Plane(vec3(0.0), vec3(0.0, 1.0, 0.0), vec3(0.5));
	Ray r = setupRay(fragCoord, 0.0);

	vec4 light = vec4(0.0);

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

	bool hitScene = traceBounce(r, plane, bg0, hit, fresnel, nml, true);

	if (hitScene) {
		float fakeBounce = 0.0;
		for (int i = 0; i < 5; i++) {
			float idx = fragCoord.y * iResolution.x * 4.0 + fragCoord.x * 4.0 + fakeBounce;
			ivec2 idxv = ivec2(mod(idx / 1024.0, 1024.0), mod(idx, 1024.0));
			if (stripes) {
				r.d = normalize(mix(reflect(r.d, nml), nml + randomVec3(idxv), (1.0-fresnel) * fract(roughness*dot(r.o, r.o)*10.0)));// + (abs(dot(r.d, nml)) * roughness) * randomVec3(5.0+r.o+r.d));
			} else {
				r.d = normalize(mix(reflect(r.d, nml), nml + randomVec3(idxv), (1.0-fresnel) * (hit.index >= 0 ? roughness : 0.05)));// + (abs(dot(r.d, nml)) * roughness) * randomVec3(5.0+r.o+r.d));
			}
			r.invD = 1.0 / r.d;
			r.o = r.o + nml * epsilon;
			if (!traceBounce(r, plane, bg0, hit, fresnel, nml, false)) {
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

	/*

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
    vec4 light = vec4(0.0);
    float error = 1.0;

    for (int k = 0; k < 8; k++) {
        if ((cameraApertureSize < 0.2 && k > 4) || (cameraApertureSize < 0.05 && k > 2)) {
            break;
        }

        Ray r = setupRay(fragCoord, float(k));

        vec3 cameraFocusVector = cameraFocusPoint - r.o;
        float focusDistance = (dot(cameraFocusVector, r.d));

        bool hitScene = traceBounce(r, plane, bg0, hit0, fresnel0, nml0, true);
        Ray primaryHit = r;

        float phi = 0.5 + 0.5 * sqrt(5.0);
        float rphi = 1.0 / phi;

        if (hitScene) {
            if (showFocalPlane) {
                r.light = primaryHit.light = r.light + 10.0 * vec3(
                    200.0 * max(0.0, 0.03 - abs(hit0.distance - focusDistance)),
                    0.0,
                    2.0 * clamp(0.1 / cameraApertureSize - abs(hit0.distance - focusDistance), 0.0, 0.1)
                );
            }

            float fakeBounce = 10.0;
            for (int j = 0; j < 8; j++) {
                if (light.a > 8.0 && error < 0.03) {
                    // light = vec4(error * 10000.0, 0.0, 0.0, 1.0);
                    return light.rgb / max(1.0, light.a);
                }
                if ((cameraApertureSize >= 0.2 && j > 2) || (cameraApertureSize >= 0.05 && j > 4)) {
                    break;
                }


                float fresnel = fresnel0;
                Hit hit = hit0;
                vec3 nml = nml0;
                for (int i = 0; i < 5; i++) {
                    if (stripes) {
                        r.d = normalize(mix(reflect(r.d, nml), nml + randomVec3(r.o+r.d+fakeBounce++), (1.0-fresnel) * mod(roughness*dot(r.o, r.o)*10.0, 1.0)));// + (abs(dot(r.d, nml)) * roughness) * randomVec3(5.0+r.o+r.d));
                    } else {
                        r.d = normalize(mix(reflect(r.d, nml), nml + randomVec3(r.o+r.d+fakeBounce++), (1.0-fresnel) * (hit.index >= 0 ? roughness : 0.05)));// + (abs(dot(r.d, nml)) * roughness) * randomVec3(5.0+r.o+r.d));
                    }
                    r.invD = 1.0 / r.d;
                    r.o = r.o + nml * epsilon;
                    if (!traceBounce(r, plane, bg0, hit, fresnel, nml, false)) {
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
    }

	return light.rgb / max(1.0, light.a);
	*/
}
