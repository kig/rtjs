
// Ray tracing setup

uniform float deviceEpsilonTrace;

uniform float cameraFocusDistance;
uniform float cameraApertureSize;

uniform bool showFocalPlane;

uniform highp usampler2D materialIndices;
uniform int materialIndicesWidth;

uniform highp sampler2D materials;
uniform int materialsWidth;

uniform highp sampler2DArray materialTextures;

const float SKY_DISTANCE = 1e6;

const vec3 texOffset = vec3(0.75, -1.25, -0.35);
const vec3 texScale = vec3(1.25);

Ray setupRay(vec2 fragCoord, float off) {
	vec2 uv = (fragCoord.xy / iResolution.xy) * 2.0 - 1.0;
	uv.x *= iResolution.x / iResolution.y;

	vec3 origin = cameraPosition;
	vec3 direction = normalize( unproject(vec3(uv, 1.0)) - origin );

	// Camera aperture simulation

	vec3 target = origin + (direction * cameraFocusDistance);

	vec3 lensOffset = applyMatrix4Rot(diskPoint(fragCoord+off) * cameraApertureSize, cameraMatrixWorld);
	origin += lensOffset;
	direction = normalize(target - origin);

	return Ray(
		origin,
		direction,
		1.0 / direction,
		vec3(1.0),
		vec3(0.0, 0.0, 0.0),
		0.0,
        1.0,
		-1
	);
}

vec4 triPlanar(vec3 nml, sampler2DArray tex, int index, vec3 p0, vec3 offset, vec3 scale) {
    vec3 p = (p0 + offset) * scale;
    // return texture(tex, toOct(p) * scale.xy, -100.0);

    vec3 blending = pow(nml, vec3(32.0));
    blending = normalize(max(blending, 0.00001)); // Force weights to sum to 1.0
    float b = (blending.x + blending.y + blending.z);
    blending /= vec3(b, b, b);

    vec4 xy = texture(tex, vec3(p.xy, index));
    vec4 xz = texture(tex, vec3(p.xz, index));
    vec4 yz = texture(tex, vec3(p.yz, index));
    return xy * blending.z + xz * blending.y + yz * blending.x;
}

vec4 sampleIndex(int index, vec3 nml, Ray r) {
    return triPlanar(nml, materialTextures, index, r.o, texOffset, texScale);
}

vec3 sampleVec3(vec3 value, vec3 nml, Ray r) {
    if (value.r < -0.5) {
        return sampleIndex(int(-value.r - 0.5), nml, r).rgb;
    } else {
        return value;
    }
}

vec3 sampleNormal(float texIndex, vec3 nml, Ray r) {
    if (texIndex < 0.0) {
        vec4 v = sampleIndex(int(-texIndex - 0.5), nml, r);
        return mix(vec3(0.0, 0.0, 1.0), v.rgb * 2.0 - 1.0, v.a);
    } else {
        return vec3(0.0, 0.0, 1.0);
    }
}

float sampleFloat(float value, vec3 nml, Ray r) {
    if (value < 0.0) {
        return sampleIndex(int(-value - 0.5), nml, r).r;
    } else {
        return value;
    }
}


Material getMaterial(Ray r, int index) {
    int materialIndex = index < 0 ? 2+index : int(readMaterialIndex(materialIndices, materialIndicesWidth, index));
    return readMaterial(materials, materialsWidth, materialIndex);
}

vec3 getTransmit(in Ray r, in int index, in Coating coating, in vec3 nml) {
    return sampleVec3(coating.color, nml, r);
}

vec3 getEmission(in Ray r, in int index, in Material material, in vec3 nml) {
    // return abs(mod(dot(r.o, r.o), 0.25)) < 0.025 ? vec3(20.0, 10.0, 3.0) : vec3(0.0);
    return sampleVec3(material.emission, nml, r);
}

vec3 getNormal(Ray r, int hitIndex, Coating coating, out vec3 rawNormal) {
    if (hitIndex < 0) {
        rawNormal = vec3(0.0, 1.0, 0.0);
        return rawNormal;
    }
    vec3 nml = triNormal(triangles, trianglesWidth, normals, normalsWidth, r.o, hitIndex);
    rawNormal = nml;

    if (coating.forwardScatter == 0.0) {
        return nml;
    }

    vec3 basisA[3];
    orthoBasis(basisA, nml);
    mat3 basis = mat3(-basisA[0], -basisA[1], basisA[2]);
    return normalize(
        basis * sampleNormal(coating.forwardScatter, nml, r)
    );
}

float fresnelP(vec3 dir, vec3 nml, float currentIOR, float nextIOR) {
    dir = normalize(dir);
    nml = normalize(nml);
    vec3 refracted = refract(dir, nml, currentIOR / nextIOR);
    float ci = dot(dir, nml);
    float ct = dot(nml, refracted);
    float Rs = (currentIOR * ci - nextIOR * ct) / (currentIOR * ci + nextIOR * ct);
    // float Rt = (currentIOR * ct - nextIOR * ci) / (currentIOR * ct + nextIOR * ci);
    Rs *= Rs;
    // Rt *= Rt;
    return Rs; //0.5 * (Rs + Rt);
}

bool getSpecular(Ray r, int hitIndex, vec3 nml, Material material, out float fresnel, out Coating coating) {
    float rn = random(r.d.xz);
    // float fresnelCoat = fresnelP(r.d, nml, r.IOR, material.coat.IOR);
    // if (rn < fresnelCoat) {
    //     coating = material.coat;
    //     fresnel = fresnelCoat;
    //     return true;
    // }
    float currentIOR = 1.0;
    float nextIOR = material.specular.IOR;
    if (dot(r.d, nml) > 0.0) {
        nextIOR = 1.0;
        currentIOR = material.specular.IOR;
        nml = -nml;
    }
    float fresnelSpecular = fresnelP(r.d, nml, currentIOR, nextIOR);
    fresnel = fresnelSpecular;
    if (rn < fresnelSpecular) {
        coating = material.specular;
        return true;
    }
    coating = material.pigment;
    return false;
}

Hit evaluateHit(in Ray r, in Plane plane) {
	Hit hit = setupHit();
	Hit hit2 = setupHit();
	intersectPlane(r, plane, hit2);
	intersectGrid(r, hit);
    if (hit2.distance < hit.distance) {
        hit = hit2;
    }
    return hit;
}

bool traceBounce(inout Ray r, in Plane plane, in vec3 bg0, out Hit hit, out float fresnel, out vec3 nml, out vec3 texNml, out Material material, out Coating coating, out vec3 transmit, in bool isPrimaryRay, inout bool specular) {
    hit = evaluateHit(r, plane);
    if (hit.distance >= SKY_DISTANCE) {
		return false;
	}
    r.light.r += float(showFocalPlane && isPrimaryRay) * 10.0 * (
        2.0 * clamp(0.1 / cameraApertureSize - abs(hit.distance - cameraFocusDistance), 0.0, 0.1) +
        200.0 * max(0.0, 0.03 - abs(hit.distance - cameraFocusDistance)) 
    );

	r.lastTested = hit.index;
	r.o += r.d * hit.distance;
	r.pathLength += hit.distance;

    material = getMaterial(r, hit.index);
	nml = getNormal(r, hit.index, material.specular, texNml);

    vec3 fog = vec3(0.0); // (1.0-exp(-hit.distance/40.0)) * bg0;
	r.light += r.transmit * (getEmission(r, hit.index, material, texNml) + fog);

    specular = getSpecular(r, hit.index, nml, material, fresnel, coating);

	transmit = getTransmit(r, hit.index, coating, texNml);

	return true;
}

vec3 skybox(in Ray r) {
    vec3 bg = clamp(dot(r.d, vec3(0.0, 1.0, 0.0)), 0.0, 1.0) * vec3(0.3, 0.5, 0.7);
    bg += vec3(12.0, 9.0, 6.0) * pow(clamp(dot(r.d, normalize(vec3(12.0, 8.0, 6.0))), 0.0, 1.0), 16.0);
    bg += vec3(10.0, 8.0, 6.0) * 4.0 * pow(clamp(dot(r.d, normalize(vec3(12.0, 8.0, 6.0))), 0.0, 1.0), 64.0);
    bg += 0.25*vec3(0.6, 0.8, 0.97) * abs(1.0+r.d.z);
    return bg;
}

vec3 trace(vec2 fragCoord, inout vec3 hitPoint, out int hitIndex, out float hitRoughness) {
	float epsilon = deviceEpsilonTrace;

	Plane plane = Plane(vec3(0.0), vec3(0.0, 1.0, 0.0), vec3(0.5));
	Ray r = setupRay(fragCoord, 0.0);

	vec4 light = vec4(0.0);

	vec3 bg0;
	bg0 = vec3(0.6+clamp(-1.0, 0.0, 1.0), 0.7, 0.8+(0.4*0.0) * abs(0.0));
	bg0 += 0.25 * vec3(10.0, 6.0, 4.0) * 4.0 * pow(clamp(dot(vec3(0.0, 1.0, 0.0), normalize(vec3(6.0, 10.0, 8.0))), 0.0, 1.0), 64.0);
	bg0 += vec3(4.0, 5.0, 7.0) * abs(1.0 - 0.0);
	bg0 *= 0.2;

	float fresnel;
	Hit hit;
	vec3 nml;
    vec3 lastNml;
    vec3 texNml;
    vec3 transmit;
    bool specular;
    Material material;
    Coating coating;

	bool hitScene = traceBounce(r, plane, bg0, hit, fresnel, nml, texNml, material, coating, transmit, true, specular);

	if (hitScene) {
        hitPoint = r.o;
        hitIndex = hit.index;
		float fakeBounce = iFrame * 1025.0;
		for (int i = 0; i < 6; i++) {
			float idx = fragCoord.y * iResolution.x * 16.0 + fragCoord.x * 16.0 + fakeBounce;
			ivec2 idxv = ivec2(mod(idx / 1024.0, 1024.0), mod(idx, 1024.0));
            float troughness = sampleFloat(coating.roughness, texNml, r);
            if (i == 0) hitRoughness = troughness;
            float randomDirFactor = (1.0-fresnel*fresnel) * troughness;

            float bounceCount = (randomDirFactor*0.1 > random(r.o.xy)) ? 2.0 : 1.0;

            vec3 retroReflectiveness = vec3(float((random(r.d.xz) < 0.1 && bounceCount > 1.0) || (mod(length(r.o), 0.1) < 0.0 && random(r.o.yz) < (0.75 - 0.75 * fresnel))));

            vec3 surfaceRay;

            float inside = sign(dot(r.d, nml));
            nml = nml * -inside;

            bool scattered = false;
            float rnd = 0.0;

            float scatterProbability = 0.0;

            if (inside > 0.0 && material.pigment.density > 0.0) {
                scatterProbability = pow(material.pigment.density, 1.0 / hit.distance);
                rnd = random(r.o.xy);
                scattered = rnd < scatterProbability;
            }

            if (scattered) {
                float scatterDistance = (1.0 - pow(random(r.o.zx), (1.0-scatterProbability))) * hit.distance;
                r.o -= r.d * (hit.distance - scatterDistance);
                r.d = normalize(mix(
                    normalize(reflect(r.d, lastNml)),
                    normalize(-r.d * (1.0 - material.pigment.forwardScatter) + randomVec3(idxv))
                , min(10.0*scatterDistance, 1.0)));
                r.transmit *= max(vec3(0.0), pow(transmit, vec3(exp(-material.pigment.density * scatterDistance)))); // pow(transmit, vec3(bounceCount));
            } else {

                if (!specular && material.pigment.density < 1.0) {
                    surfaceRay = refract(r.d, nml, inside < 0.0 ? (1.0 / material.specular.IOR) : (material.specular.IOR / 1.0));
                    nml = -nml;
                    r.IOR = inside < 0.0 ? material.specular.IOR : 1.0;
                    //r.transmit *= max(vec3(0.0), transmit);
                } else {
                    surfaceRay = reflect(r.d, nml);
                    r.transmit *= max(vec3(0.0), pow(transmit, vec3(bounceCount)));
                }

                r.d = normalize(mix(mix(
                    surfaceRay, 
                    nml + randomInUnitSphere(idxv),
                    randomDirFactor
                ), -r.d + 0.1*randomVec3(idxv), retroReflectiveness));
                r.o = r.o + nml * epsilon;
                lastNml = nml;
            }
            r.invD = 1.0 / r.d;
            if (length(transmit) < 0.01) {
                break;
            }
			if (!traceBounce(r, plane, bg0, hit, fresnel, nml, texNml, material, coating, transmit, false, specular)) {
				light += vec4(max(vec3(0.0), r.light + r.transmit * skybox(r)), 1.0);
				break;
			}
            fakeBounce++;
		}
	} else {
        light += vec4(r.light + (float(!costVis) * r.transmit) * skybox(r), 1.0);
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
