
// Ray tracing setup

uniform float deviceEpsilonTrace;

uniform float cameraFocusDistance;
uniform float cameraApertureSize;

uniform float roughness;

uniform bool showFocalPlane;
uniform bool stripes;

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

    vec3 blending = pow(nml, vec3(8.0));
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
    int materialIndex = 0;
    if (index >= 0) {
        materialIndex = int(readMaterialIndex(materialIndices, materialIndicesWidth, index));
    }
    return readMaterial(materials, materialsWidth, materialIndex);
}

vec3 getTransmit(in Ray r, in int index, in Coating coating, in vec3 nml) {
	if (index < 0) {
		return vec3(0.05);
	} else {
        return sampleVec3(coating.color, nml, r);
	}
}

vec3 getEmission(in Ray r, in int index, in Material material, in vec3 nml) {
    if (index < 0) {
        return mix(vec3(0.0), vec3(10.0), float(length(r.o) > 1.2) * clamp(1.0 - (length(r.o) - 1.2) / 0.4, 0.0, 1.0) );
    }
    return sampleVec3(material.emission, nml, r);
}

vec3 getNormal(Ray r, int hitIndex, Coating coating, out vec3 rawNormal) {
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

bool getSpecular(Ray r, int hitIndex, Material material) {
    return false;
}

bool traceBounce(inout Ray r, in Plane plane, in vec3 bg0, out Hit hit, out float fresnel, out vec3 nml, out vec3 texNml, out Material material, out vec3 transmit, in bool isPrimaryRay, inout bool specular) {
	hit = setupHit();
	Hit hit2 = setupHit();
	// intersectSphere(r, vec3(0.0, 2.5, 0.0), 0.5, hit);
	intersectGrid(r, hit);
	intersectPlane(r, plane, hit2);
	if (hit2.distance < hit.distance) {
		hit = hit2;
	}
	if (hit.distance >= SKY_DISTANCE) {
		return false;
	}
    // r.light.r += float(showFocalPlane && isPrimaryRay) * 10.0 * (
    //     2.0 * clamp(0.1 / cameraApertureSize - abs(hit.distance - cameraFocusDistance), 0.0, 0.1) +
    //     200.0 * max(0.0, 0.03 - abs(hit.distance - cameraFocusDistance)) 
    // );

	r.lastTested = hit.index;
	r.o += r.d * hit.distance;
	r.pathLength += hit.distance;

    material = getMaterial(r, hit.index);
    specular = getSpecular(r, hit.index, material);
    Coating coating = material.pigment;
    if (specular) {
        coating = material.specular;
    }
	nml = hit.index >= 0 
        ? getNormal(r, hit.index, coating, texNml)
        : plane.normal;

    vec3 fog = (1.0-exp(-hit.distance/40.0)) * bg0;
	r.light += r.transmit * (getEmission(r, hit.index, material, texNml) + fog);

	fresnel = pow(1.0 - abs(dot(r.d, nml)), 5.0);
	transmit = getTransmit(r, hit.index, coating, texNml);

	return true;
}

vec3 skybox(in Ray r) {
    vec3 bg = vec3(0.6+clamp(-r.d.y, 0.0, 1.0), 0.7, 0.8+(0.4*r.d.x) * abs(r.d.z));
    bg += vec3(10.0, 9.0, 7.0) * pow(clamp(dot(r.d, normalize(vec3(6.0, 10.0, 8.0))), 0.0, 1.0), 16.0);
    bg += vec3(10.0, 8.0, 6.0) * 4.0 * pow(clamp(dot(r.d, normalize(vec3(6.0, 10.0, 8.0))), 0.0, 1.0), 64.0);
    bg += 0.5 * 5.5*vec3(0.7, 0.8, 1.0) * abs(1.0+r.d.y);
    return bg;
}

vec3 trace(vec2 fragCoord) {
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
    vec3 texNml;
    vec3 transmit;
    bool specular;
    Material material;

	bool hitScene = traceBounce(r, plane, bg0, hit, fresnel, nml, texNml, material, transmit, true, specular);

	if (hitScene) {
		float fakeBounce = iFrame * 7.0 * 1025.0;
        bool sss = false;
		for (int i = 0; i < 4; i++) {
			float idx = fragCoord.y * iResolution.x * 16.0 + fragCoord.x * 16.0 + fakeBounce;
			ivec2 idxv = ivec2(mod(idx / 1024.0, 1024.0), mod(idx, 1024.0));
            float troughness = sampleFloat(specular ? material.specular.roughness : material.pigment.roughness, texNml, r);
            float randomDirFactor = (1.0-fresnel*fresnel)*(hit.index >= 0 ? troughness : fract(0.5*dot(r.o, r.o)*10.0));

            float bounceCount = 1.0; // (randomDirFactor*0.1 > random(r.o.xy)) ? 2.0 : 1.0;

            // vec3 retroReflectiveness = hit.index >= 0
            //     ? vec3(float((random(r.d.xz) < 0.1 && bounceCount > 1.0) || (mod(length(r.o), 0.1) < 0.0 && random(r.o.yz) < (0.75 - 0.75 * fresnel))))
            //     : vec3(0.0);

            vec3 surfaceRay;

            float inside = sign(dot(r.d, nml));
            nml = nml * -inside;

            // float scatterDistance = pow(random(r.d.xz), 0.5);

            // if (hit.index >= 0 && inside > 0.0 && hit.distance > scatterDistance && !sss) {
            //     r.o -= r.d * (hit.distance - scatterDistance);
            //     r.d = normalize(-r.d * 0.8 + randomVec3(idxv));
            //     bounceCount = 0.5;
            //     sss = true;
            // } else {
                sss = false;

                // if (!specular && hit.index >= 0) {
                //     surfaceRay = refract(r.d, nml, inside < 0.0 ? 1.0/1.84 : 1.84/1.0);
                //     nml = -nml;
                // } else {
                    surfaceRay = reflect(r.d, nml);
                // }

                r.d = normalize(mix(
                    surfaceRay, 
                    nml + randomVec3(idxv),
                    randomDirFactor
                )); //, -r.d + randomVec3(idxv), retroReflectiveness)); // + (abs(dot(r.d, nml)) * roughness) * randomVec3(5.0+r.o+r.d));
                r.o = r.o + nml * epsilon;
            // }
            r.transmit *= pow(transmit, vec3(bounceCount));
            r.invD = 1.0 / r.d;
			if (!traceBounce(r, plane, bg0, hit, fresnel, nml, texNml, material, transmit, false, specular)) {
				light += vec4(r.light + (float(!costVis) * r.transmit) * skybox(r), 1.0);
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
