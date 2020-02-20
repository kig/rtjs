const mobile = /mobile/i.test(navigator.userAgent);
var dpr = mobile ? 1 : (window.devicePixelRatio || 1);
var resAdjust = 1;
if (/Mac OS X/.test(navigator.userAgent)) {
    if (window.screen.width !== 1280 && dpr === 2) {
        resAdjust = 1280 / window.screen.width;
    }
}
dpr *= resAdjust;

class Coating {

    constructor(options) {
        this.color = [0,0,0];
        this.roughness = 0;
        this.IOR = 1;
        this.retroReflectivity = 0.0;
        this.anisotropy = 0;
        this.anisotropyAngle = 0;
        this.density = 1;
        this.forwardScatter = 0.5;
        this.thickness = 0;
        this.transmission = 0;
        this.subsurface = 0;

        if (options) {
            for (let name in options) {
                this[name] = options[name];
            }
        }
    }

    toArray() {
        return this.color.concat([
            this.roughness, this.IOR, this.retroReflectivity,
            this.anisotropy, this.anisotropyAngle,
            this.density, this.forwardScatter,
            this.thickness,
            this.transmission,
            this.subsurface
        ]);
    }

}

class Material {

    constructor(options) {
        this.specular = new Coating({color: [1,1,1], IOR: 1.5});
        this.coat = new Coating();
        this.volume = new Coating({roughness: 1, IOR: 10});
        this.thinWall = 0;
        this.wallThickness = 0;
        this.emission = [0,0,0];
        if (options) {
            for (let name in options) {
                this[name] = options[name];
            }
        }
    }

    toArray() {
        return this.specular.toArray().concat(
            this.coat.toArray(), 
            this.volume.toArray(),
            [this.thinWall, this.wallThickness],
            this.emission
        );
    }

}

Material.makeRandom = function() {
    const m = new Material();
    m.specular.color = vecToArray(addS(mulS(randomVec3Positive(), 0.05), 0.95));
    m.coat.color = vecToArray(randomVec3Positive());
    m.volume.color = vecToArray(addS(mulS(randomVec3Positive(), 0.8), 0.1));
    // if (random() < 0.01) {
    //     m.emission = vecToArray(addS(mulS(randomVec3Positive(), 40.15), 2.05));
    // }
    m.volume.roughness = random();
    m.specular.roughness = random();
    m.specular.IOR = 1 + random();
    m.volume.density = 1;
    m.coat.IOR = 1 + random();
    return m.toArray();
};

Material.glassy = (function() {
    const m = new Material();
    m.specular.color = vecToArray(vec3(1.0));
    m.coat.color = vecToArray(vec3(0.85, 0.51, 0.35));
    m.volume.color = vecToArray(vec3(0.9, 0.71, 0.5));
    m.volume.roughness = 0.0;
    m.specular.roughness = 0.05;
    m.specular.IOR = 1.5;
    m.volume.density = 0.05;
    m.coat.IOR = 1;
    return m.toArray();
})();

Material.steel = (function() {
    var specular = [
        1,1,1, // color
        -3,    // roughness
        1.5,   // IOR
        0.05,  // retroreflectiveness
        0, 0,  // anisotropy
        -2,     // density
        -4,     // forwardScatter
        0,     // thickness
        0,     // transmission
        0      // subsurface
    ];
    var coat = [
        0,1,0, // color
        0,   // roughness
        1.5,    // IOR
        0,  // retroreflectiveness
        0, 0,  // anisotropy
        0,     // density
        0,     // forwardScatter
        0,     // thickness
        0,     // transmission
        0      // subsurface
    ];
    var volume = [
        -1,0.5,0.5, // color
        -3,   // roughness
        10,    // IOR
        0.05,  // retroreflectiveness
        0, 0,  // anisotropy
        1,     // density
        -4,     // forwardScatter
        0,     // thickness
        0,     // transmission
        0      // subsurface
    ];
    var emission = [0,0,0];
    return specular.concat(coat).concat(volume).concat([0,0]).concat(emission);
})();


class WebGLTracer2 {
    constructor(vg, traceGLSL, blueNoiseTexture, diffuseTexture, metallicTexture, roughnessTexture, normalTexture) {
        var canvas = document.createElement( 'canvas' );
        var context = canvas.getContext( 'webgl2' );

        this.diffuseTexture = diffuseTexture;
        this.metallicTexture = metallicTexture;
        this.roughnessTexture = roughnessTexture;
        this.normalTexture = normalTexture;

        this.voxelGrid = vg;

        console.time('Serialize VoxelGrid');
        const arrays = vg.serialize();
        this.triangleCount = arrays.triangles.length;
        console.timeEnd('Serialize VoxelGrid');

        this.stats = {
            element: document.createElement('div'),
            fields: {},
            log: function(name, msg) {
                if (!this.fields[name]) {
                    this.fields[name] = document.createElement('div');
                    this.fields[name].titleEl = document.createElement('span');
                    this.fields[name].valueEl = document.createElement('span');
                    this.fields[name].titleEl.textContent = name;
                    this.fields[name].append(this.fields[name].titleEl);
                    this.fields[name].append(this.fields[name].valueEl);
                    this.element.append(this.fields[name]);
                }
                this.fields[name].valueEl.textContent = msg;
            }
        };
        this.stats.element.className = 'stats';

        document.body.querySelector('#controls').appendChild(this.stats.element);

        const toBytes = (arr) => {
            const u8 = new Uint8Array(arr.length * 32);
            for (let i=0; i<u8.length; i++) {
                u8[i] = (arr[i >> 5] >> (i & 31)) & 1;
            }
            return u8;
        };

        

        var triCount = [];
        for (let i = 0, j = 0; i < this.triangleCount; i += 9, j++) {
            triCount.push(1);
        }
        var materials = [];
        Material.makeRandom().forEach(k => materials.push(k));
        Material.steel.forEach(k => materials.push(k));
        for (let i = 0; i < triCount.length; i++) {
            Material.makeRandom().forEach(k => materials.push(k));
        }



        this.textures = {
            voxelIndex: this.createTexture(new Int16Array(arrays.voxelIndex.buffer), THREE.RedIntegerFormat, THREE.ShortType),
            triIndices: this.createTexture(arrays.triIndices, THREE.RedIntegerFormat, THREE.UnsignedIntType),
            triangles: this.createTexture(arrays.triangles, THREE.RedFormat, THREE.FloatType),
            normals: this.createTexture(arrays.normals, THREE.RedFormat, THREE.FloatType),
            materialIndices: this.createTexture(new Uint32Array(triCount), THREE.RedIntegerFormat, THREE.UnsignedIntType),
            materials: this.createTexture(new Float32Array(materials), THREE.RedFormat, THREE.FloatType)
        };

        var materialImages = [
            this.normalTexture.image,
            this.roughnessTexture.image,
            this.metallicTexture.image,
            this.diffuseTexture.image
        ];

        var tcanvas = document.createElement('canvas');
        var tctx = tcanvas.getContext('2d');
        var texWidth = materialImages[0].width;
        var texHeight = materialImages[0].height;
        var texDepth = materialImages.length;
        tcanvas.width = texWidth;
        tcanvas.height = texHeight * texDepth;
        tctx.translate(0, tcanvas.height);
        tctx.scale(1, -1);
        materialImages.forEach((img, i) => {
            tctx.drawImage(img, 0, i * texHeight);
        });
        var imageData = tctx.getImageData(0, 0, tcanvas.width, tcanvas.height);

        tcanvas.setAttribute('style', 'right: 0; z-index: 10; left: auto; width: auto !important;');
        // document.body.appendChild(tcanvas);

        this.materialTextures = new THREE.DataTexture3D(
            new Uint8Array(imageData.data.buffer), texWidth, texHeight, texDepth
        );
        this.materialTextures.wrapS = THREE.RepeatWrapping;
        this.materialTextures.wrapT = THREE.RepeatWrapping;
        this.materialTextures.wrapR = THREE.RepeatWrapping;
        this.materialTextures.needsUpdate = true;

        this.renderer = new THREE.WebGLRenderer({ canvas, context });
        this.renderer.setSize(window.innerWidth*dpr, window.innerHeight*dpr);
        this.renderer.setClearColor(0x00ffff);
        this.renderer.clear();

        this.renderTarget = new THREE.WebGLRenderTarget(window.innerWidth*dpr, window.innerHeight*dpr, {
            format: THREE.RGBAFormat,
            type: THREE.FloatType
        });
        this.accumRenderTargetA = new THREE.WebGLRenderTarget(window.innerWidth*dpr, window.innerHeight*dpr, {
            format: THREE.RGBAFormat,
            type: THREE.FloatType
        });
        this.accumRenderTargetB = new THREE.WebGLRenderTarget(window.innerWidth*dpr, window.innerHeight*dpr, {
            format: THREE.RGBAFormat,
            type: THREE.FloatType
        });
        this.varianceRenderTarget = new THREE.WebGLRenderTarget(window.innerWidth*dpr, window.innerHeight*dpr, {
            format: THREE.RGBAFormat,
            type: THREE.FloatType
        });

        var camera = new THREE.PerspectiveCamera(55, 1, 0.1, 100);
        camera.setFocalLength(50);
        camera.target = new THREE.Vector3(0, 1.1, 0);
        this.focusPoint = new THREE.Vector3(-0.9, 1.3, 0.3);
        this.focusDistance = this.focusPoint.distanceTo(camera.position);
		camera.apertureSize = Math.pow(1.33, -11);
        camera.previousMatrix = new THREE.Matrix4();
        camera.uvToWorld = new THREE.Matrix4();
        camera.worldToUV = new THREE.Matrix4();
        camera.lookAt(camera.target);
        camera.updateProjectionMatrix();
        camera.updateMatrixWorld();

        this.material = new THREE.RawShaderMaterial({
            uniforms: {
                iTime: { value: 1.0 },
                iFrame: { value: 0 },
                rayBudget: { value: 1 },

                triangles: { value: this.textures.triangles },
                trianglesWidth: { value: this.textures.triangles.image.width },
                normals: { value: this.textures.normals },
                normalsWidth: { value: this.textures.normals.image.width },
                voxelIndex: { value: this.textures.voxelIndex },
                voxelIndexWidth: { value: this.textures.voxelIndex.image.width },
                triIndices: { value: this.textures.triIndices },
                triIndicesWidth: { value: this.textures.triIndices.image.width },
                materialIndices: { value: this.textures.materialIndices },
                materialIndicesWidth: { value: this.textures.materialIndices.image.width },
                materials: { value: this.textures.materials },
                materialsWidth: { value: this.textures.materials.image.width },

                materialTextures: {
                    value: this.materialTextures
                },

                varianceTexture: { value: this.varianceRenderTarget.texture },
                tex: { value: this.renderTarget.texture },

                vgOrigin: { value: new THREE.Vector3().copy(vg._origin) },
                vgScale: { value: vg._scale },
                vgSize: { value: vg._size },

                blueNoise: { value: blueNoiseTexture },

                iResolution: { value: [canvas.width, canvas.height] },
                cameraFocusDistance: { value: this.focusDistance },
                cameraApertureSize: { value: camera.apertureSize },
                uvToWorld: { value: camera.uvToWorld },
                cameraMatrixWorld: { value: camera.matrixWorld },
                previousCameraMatrixWorld: { value: new THREE.Matrix4() },
                previousCameraProjectionMatrix: { value: new THREE.Matrix4() },
                previousCameraPosition: { value: new THREE.Vector3() },
                deviceEpsilon: {value: mobile ? 0.01 : 0.001},
                deviceEpsilonTrace: {value: mobile ? 0.05 : 0.01},

                useTemporalReprojection: { value: true },
                temporalReprojectionShowRejection: { value: false },
                temporalReprojectionVarianceCutoff: { value: 2.0 },
                temporalReprojectionWeight: { value: 20.0 },

                firstFrameSampleBoost: {value: 1.0},

                costVis: {value: false},
                aaSize: {value: 4.0},
                showFocalPlane: { value: false },
                showBoost: { value: false }
            },
            vertexShader: `#version 300 es

            precision highp float;
            precision highp int;

            in vec3 position;

            void main() {
                gl_Position = vec4(position.xyz, 1.0);
            }
            `,
            fragmentShader: `#version 300 es

            precision highp float;
            precision highp int;

            uniform vec3 cameraPosition;
            uniform vec3 previousCameraPosition;
            uniform mat4 previousCameraMatrixWorld;
            uniform mat4 previousCameraProjectionMatrix;

            uniform sampler2D varianceTexture;
            uniform sampler2D tex;

            uniform bool useTemporalReprojection;
            uniform bool temporalReprojectionShowRejection;
            uniform float temporalReprojectionVarianceCutoff;
            uniform float temporalReprojectionWeight;

            uniform float firstFrameSampleBoost;

            ${traceGLSL}

            uniform float aaSize;
            uniform float rayBudget;
            uniform bool showBoost;

            out vec4 FragColor;

            void main() {

                float samples = 1.0;
                vec4 sum = vec4(0.0);

                vec2 centerToPixel = vec2(iResolution.x / iResolution.y, 1.0) * (gl_FragCoord.xy / iResolution.xy - 0.5);
                vec2 centerToPixel8x8 = vec2(iResolution.x / iResolution.y, 1.0) * (8.0*floor(gl_FragCoord.xy / 8.0) / iResolution.xy - 0.5);
                float distanceToCenterPx = length(centerToPixel);
                float distanceToCenter = length(centerToPixel8x8);
                float boostRadius = max(0.2, rayBudget * 0.4);
                float fullBoostRadius = max(0.025, rayBudget * 0.1);
                float boostExponent = 8.0 / max(0.125, rayBudget);
                float fullBoostSamples = rayBudget * 10.0;
                float boostFactor = (1.0 - clamp((distanceToCenter - fullBoostRadius) / boostRadius, 0.0, 1.0));
                boostFactor *= 1.0 - clamp(iFrame / 10.0, 0.0, 1.0);

                vec4 varianceMetrics = texelFetch(varianceTexture, ivec2(gl_FragCoord.xy), 0);
                float errorLum = varianceMetrics.x;
                float totalVariance = varianceMetrics.y;
                bool converged = varianceMetrics.z > 0.0;
                bool convergedVariance = varianceMetrics.w > 0.0;

                if (iFrame > 1.0) {

                    if (( ((iFrame > 2.0 && iFrame < 10.0) || (iFrame > 30.0 && errorLum < 0.0001 && totalVariance < 0.01)) && (converged || convergedVariance || (iFrame >= 2.0 && (errorLum < 0.01 || totalVariance < 0.01))) )) {
                        FragColor = vec4(0.0, 0.0, 0.0, 0.0);
                        return;
                    }

                    boostFactor += totalVariance/4.0;
                    boostFactor = clamp(boostFactor, 0.0, 1.0);
                    float boost = pow(boostFactor, boostExponent) * fullBoostSamples;
                    float sampleCountJitter = 0.0;
                    samples = boost + sampleCountJitter;

                    if (iFrame > 10.0) {
                        samples += totalVariance + fract(random(vec2(IR_2+iTime, IR_1*0.1+iTime)));
                    }
                    
                } else {
                    boostFactor = clamp(boostFactor, 0.0, 1.0);
                    float boost = pow(boostFactor, boostExponent) * fullBoostSamples;
                    float sampleCountJitter = 0.0;
                    samples = max(1.0, boost + sampleCountJitter);
                }
                if (iFrame == 1.0) samples *= firstFrameSampleBoost;

                if (showBoost) {
                    sum.g += samples;
                    if (distanceToCenterPx < fullBoostRadius) {
                        if (distanceToCenterPx < 0.002 || (fullBoostRadius - distanceToCenterPx) < 0.002) {
                            sum.r += samples;
                            sum.gb *= 0.0;
                        }
                    }
                }

                vec3 hitPoint;
                int hitIndex;
                float hitRoughness;
                bool fetched = false;
                bool applied = false;
                for (float i = 0.0; i < samples; i++) {
                    hitPoint = vec3(1.0e6+1.0);
                    float y = mod(i, 2.0);
                    float x = i - y * 2.0;
                    float ry = mod(y + iFrame / aaSize, aaSize);
                    float rx = mod(x + iFrame - aaSize * ry, aaSize);
                    vec3 c = trace(gl_FragCoord.xy + (vec2(rx,ry) + vec2(random(i+0.5*vec2(rx,ry)), random(i+0.5+0.5*vec2(rx,ry)))) / aaSize, hitPoint, hitIndex, hitRoughness);
                    if (useTemporalReprojection && !fetched && hitPoint.x < 1.0e6 && iFrame < 1.0) {
                        vec4 uv = previousCameraProjectionMatrix * inverse(previousCameraMatrixWorld) * vec4(hitPoint, 1.0);
                        uv /= uv.w;
                        uv.x /= (iResolution.x / iResolution.y);
                        uv = (uv + 1.0) / 2.0;
                        if (all(greaterThanEqual(uv.xy, vec2(0.0))) && all(lessThanEqual(uv.xy, vec2(1.0)))) {
                            fetched = true;
                            vec3 direction = normalize(hitPoint - previousCameraPosition);
                            Ray r = Ray(
                                previousCameraPosition,
                                direction,
                                1.0 / direction,
                                vec3(1.0),
                                vec3(0.0, 0.0, 0.0),
                                0.0,
                                1.0,
                                -1                        
                            );
                            Hit hit = evaluateHit(r, Plane(vec3(0.0), vec3(0.0, 1.0, 0.0), vec3(0.5)));
                            if (hit.distance < SKY_DISTANCE && hit.index == hitIndex) {
                                vec4 previousVariance = texture(varianceTexture, uv.xy, 0.0);
                                if (previousVariance.x < temporalReprojectionVarianceCutoff) {
                                    vec4 start = texture(tex, uv.xy, 0.0);
                                    sum += vec4(start.rgb / max(1.0, start.a) * min(start.a, temporalReprojectionWeight), min(start.a, temporalReprojectionWeight));
                                    applied = true;
                                } else if (temporalReprojectionShowRejection) {
                                    sum += vec4(0.0, 400.0, 0.0, 100.0);
                                }
                            } else if (temporalReprojectionShowRejection) {
                                sum += vec4(400.0, 0.0, 0.0, 100.0);
                            }
                        }
                    }
                    if (applied && c.g > (sum.g/sum.a)/0.75) c *= 0.75;
                    sum += vec4(c, 1.0);
                }
                if (applied) sum = vec4((sum.rgb / sum.a) * (samples+(sum.a-samples)*0.8), (samples+(sum.a-samples)*0.8));
                FragColor = sum;
            }
            `,
            depthTest: false,
            depthWrite: false
        });
        this.mesh = new THREE.Mesh(
            new THREE.PlaneGeometry(2,2,1),
            this.material
        );
        this.mesh.frustumCulled = false;
        this.mesh.position.z = -0.5;
        this.mesh.rotation.y = Math.PI;

        this.varianceMaterial = new THREE.RawShaderMaterial({
            uniforms: {
                tex: { value: this.accumRenderTargetB.texture },
                previousTex: { value: this.renderTarget.texture },
                iFrame: { value: this.frame }
            },
            vertexShader: this.material.vertexShader,
            fragmentShader: `#version 300 es

            precision highp float;
            precision highp int;
            
            uniform sampler2D tex;
            uniform sampler2D previousTex;

            uniform float iFrame;

            out vec4 FragColor;    

            float rgbToPerceivedLuminance(vec3 c) {
                return clamp(c.r * 0.2126 + c.g * 0.7152 + c.b * 0.0722, 0.0, 1.0);
            }

            vec3 bilinear(vec3 tl, vec3 tr, vec3 bl, vec3 br, vec2 uv) {
                return mix(mix(tl, tr, uv.x), mix(bl, br, uv.x), uv.y);
            }

            vec3 toGamma(vec4 c) {
                return clamp((1.0 - exp(-0.5 * c.rgb / c.a)), vec3(0.0), vec3(1.0));
            }

            void main() {

                // if (iFrame == 0.0) {
                //     FragColor = vec4(1.0, 1.0, 0.0, 0.0);
                //     return;
                // }

                vec4 prev = texelFetch(previousTex, ivec2(gl_FragCoord.xy), 0);
                vec4 current = texelFetch(tex, ivec2(gl_FragCoord.xy), 0);
                vec3 error = toGamma(current) - toGamma(prev);
                float errorLum = rgbToPerceivedLuminance(abs(error));

                vec4 pixels[16];
                bool converged = true;
                vec4 sum = vec4(0.0);
                for (int y = 0, i = 0; y < 4; y++)
                for (int x = 0; x < 4; x++, i++) {
                    vec4 px = texelFetch(tex, ivec2(4.0 * floor(gl_FragCoord.xy/4.0)) + ivec2(x, y), 0);
                    pixels[i] = vec4(toGamma(px), px.a);
                    converged = converged && (px.a >= 500.0);
                    sum += px;
                }
                bool convergedVariance = true;
                vec3 avg = sum.rgb / sum.a;
                float totalVariance = 0.0;
                vec3 tl = pixels[0].rgb;
                vec3 tr = pixels[3].rgb;
                vec3 bl = pixels[12].rgb;
                vec3 br = pixels[15].rgb;
                for (int y = 0, i = 0; y < 4; y++)
                for (int x = 0; x < 4; x++, i++) {
                    vec4 px = pixels[i];
                    vec3 v = px.rgb;
                    vec3 g = bilinear(tl, tr, bl, br, vec2(float(x)/3.0, float(y)/3.0));
                    vec3 d = v - g;
                    float variance = rgbToPerceivedLuminance(abs(v - g));
                    convergedVariance = convergedVariance && (px.a >= 16.0 && variance < (0.5 / 256.0));
                    totalVariance += variance;
                }

                FragColor = vec4(errorLum, totalVariance, float(converged), float(convergedVariance));
            }
            `,
            depthTest: false,
            depthWrite: false
        });
        this.varianceMesh = new THREE.Mesh(
            new THREE.PlaneGeometry(2,2,1),
            this.varianceMaterial
        );
        this.varianceMesh.frustumCulled = false;
        this.varianceMesh.position.z = -0.5;
        this.varianceMesh.rotation.y = Math.PI;

        this.accumMaterial = new THREE.RawShaderMaterial({
            uniforms: {
                tex: { value: this.renderTarget.texture },
                accumTex: { value: this.accumRenderTargetA.texture },
                iFrame: { value: 0 }
            },
            vertexShader: this.material.vertexShader,
            fragmentShader: `#version 300 es

            precision highp float;
            precision highp int;
            
            uniform sampler2D tex;
            uniform sampler2D accumTex;
            
            uniform float iFrame;

            out vec4 FragColor;    

            void main() {
                vec4 src = texelFetch(tex, ivec2(gl_FragCoord.xy), 0);
                vec4 dst = texelFetch(accumTex, ivec2(gl_FragCoord.xy), 0);
                FragColor = iFrame > 0.0 ? dst + src : src;
            }
            `,
            depthTest: false,
            depthWrite: false
        });
        this.accumMesh = new THREE.Mesh(
            new THREE.PlaneGeometry(2,2,1),
            this.accumMaterial
        );
        this.accumMesh.frustumCulled = false;
        this.accumMesh.position.z = -0.5;
        this.accumMesh.rotation.y = Math.PI;


        this.blurMaterial = new THREE.RawShaderMaterial({
            uniforms: {
                tex: { value: this.accumRenderTargetB.texture },
                sigma: {value: this.frame+1 },
                direction: { value: 0 }
            },
            vertexShader: this.material.vertexShader,
            fragmentShader: `#version 300 es

            precision highp float;
            precision highp int;
            
            uniform sampler2D tex;
            uniform float sigma;
            uniform int direction;

            out vec4 FragColor;    

            float normpdf(in float x, in float sigma) {
                return 0.39894 * exp(-0.5 * x * x / (sigma * sigma)) / sigma;
            }

            void main() {
                const int radius = 30;
                
                vec4 accum = normpdf(0.0, sigma) * texelFetch(tex, ivec2(gl_FragCoord.xy), 0);

                for (int i = 1; i <= radius; i++) {
                    ivec2 offset = ivec2(i * direction, i * (1 - direction));
                    accum += normpdf(float(i), sigma) * (
                        texelFetch(tex, ivec2(gl_FragCoord.xy) - offset, 0) +
                        texelFetch(tex, ivec2(gl_FragCoord.xy) + offset, 0)
                    );
                }

                FragColor = accum;
            }
            `,
            depthTest: false,
            depthWrite: false
        });
        this.blurMesh = new THREE.Mesh(
            new THREE.PlaneGeometry(2,2,1),
            this.blurMaterial
        );
        this.blurMesh.frustumCulled = false;
        this.blurMesh.position.z = -0.5;
        this.blurMesh.rotation.y = Math.PI;

        this.bloomMaterial = new THREE.RawShaderMaterial({
            uniforms: {
                tex: { value: this.accumRenderTargetB.texture },
                sigma: {value: 7.5 },
                startBrightness: {value: 4},
                softKnee: {value: 0},
                intensity: {value: 1}
            },
            vertexShader: this.material.vertexShader,
            fragmentShader: `#version 300 es

            precision highp float;
            precision highp int;
            
            uniform sampler2D tex;
            uniform float sigma;

            uniform float startBrightness;
            uniform float softKnee;
            uniform float intensity;

            out vec4 FragColor;    

            float normpdf(in float x, in float sigma) {
                return 0.39894 * exp(-0.5 * x * x / (sigma * sigma)) / sigma;
            }

            void main() {
                const float radius = 1.0;
                
                vec4 accum = texelFetch(tex, ivec2(gl_FragCoord.xy), 0);
                accum.rgb /= accum.a;

                for (float y = -radius; y <= radius; y++) {
                    for (float x = -radius; x <= radius; x++) {
                        if (x == 0.0 && y == 0.0) {
                            continue;
                        }
                        float d = pow(length(vec2(x,y)/2.0), 6.0);
                        vec4 c = texelFetch(tex, ivec2(gl_FragCoord.xy + vec2(x,y)), 0);
                        c.rgb /= c.a;
                        if ((c.r+c.g+c.b)/3.0 > startBrightness) {
                            accum.rgb += normpdf(d, sigma) * max(vec3(0.0), c.rgb - startBrightness) * intensity;
                        }
                    }
                }

                accum.rgb *= accum.a;

                FragColor = accum;
            }
            `,
            depthTest: false,
            depthWrite: false
        });
        this.bloomMesh = new THREE.Mesh(
            new THREE.PlaneGeometry(2,2,1),
            this.bloomMaterial
        );
        this.bloomMesh.frustumCulled = false;
        this.bloomMesh.position.z = -0.5;
        this.bloomMesh.rotation.y = Math.PI;

        this.blitMaterial = new THREE.RawShaderMaterial({
            uniforms: {
                tex: { value: this.renderTarget.texture },
                varianceTexture: { value: this.varianceRenderTarget.texture },
                showConverged: { value: false },
                iFrame: { value: 0 },
                showSampleCount: { value: false }
            },
            vertexShader: this.material.vertexShader,
            fragmentShader: `#version 300 es

            precision highp float;
            precision highp int;
            
            uniform sampler2D tex;
            uniform sampler2D varianceTexture;

            uniform float iFrame;

            uniform bool showConverged;
            uniform bool showSampleCount;

            out vec4 FragColor;

            vec3 toGamma(vec4 c) {
                return clamp((1.0 - exp(-0.5 * c.rgb / c.a)), vec3(0.0), vec3(1.0));
            }
            
            void main() {
                vec4 src = texelFetch(tex, ivec2(gl_FragCoord.xy), 0);
                FragColor = vec4(toGamma(src), 1.0);
                vec4 varianceMetrics = texelFetch(varianceTexture, ivec2(gl_FragCoord.xy), 0);
                float errorLum = varianceMetrics.x;
                float totalVariance = varianceMetrics.y;
                bool converged = varianceMetrics.z > 0.0;
                bool convergedVariance = varianceMetrics.w > 0.0;

                if (showConverged) {
                    if (( ((errorLum < 0.0001)) && (converged || convergedVariance || (iFrame >= 2.0 && (errorLum < 0.01 || totalVariance < 0.01))) )) {
                        FragColor.r -= 0.5;
                        FragColor.b += 0.5;
                    } else {
                        FragColor.r += 0.5;
                        FragColor.b -= 0.5;
                    }
                }
                if (showSampleCount) {
                    float f = mod(src.a/15000.0, 1.0);
                    FragColor = vec4(
                        f < 0.25 
                        ? mix(vec3(0.1, 0.2, 0.5), vec3(0.5, 0.8, 0.2), vec3(sqrt(f*4.0)))
                        : mix(vec3(0.5, 0.8, 0.2), vec3(1.0, 0.4, 0.0), vec3((f-0.25)/0.75*4.0)), 
                        1.0);
                }
            }
            `,
            depthTest: false,
            depthWrite: false
        });
        this.blitMesh = new THREE.Mesh(
            new THREE.PlaneGeometry(2,2,1),
            this.blitMaterial
        );
        this.blitMesh.frustumCulled = false;
        this.blitMesh.position.z = -0.5;
        this.blitMesh.rotation.y = Math.PI;

        this.scene = new THREE.Scene();
        this.camera = camera;

        this.scene.add(this.camera);
        this.scene.add(this.mesh);

        this.frame = 0;

        this.controls = new CameraControls(camera, canvas);
        camera.positionOffset.y = 0.1;
        this.startTime = Date.now();

        this.rayBudget = 1;
        this.startRayBudget = 1;
        this.frameTime = 0;
        this.frameStartTime = Date.now();

        window.onresize = () => {
            this.controls.changed = true;
        };
    }

    useRandomMaterials() {
        var triCount = [];
        for (let i = 0, j = 0; i < this.triangleCount; i += 9, j++) {
            triCount.push(j+1);
        }
        var materials = [];
        Material.makeRandom().forEach(k => materials.push(k));
        for (let i = 0; i < triCount.length; i++) {
            Material.makeRandom().forEach(k => materials.push(k));
        }
        this.setMaterials(triCount, materials);
    }

    useSteelMaterial() {
        var triCount = [];
        for (let i = 0, j = 0; i < this.triangleCount; i += 9, j++) {
            triCount.push(1);
        }
        var materials = [];
        Material.makeRandom().forEach(k => materials.push(k));
        Material.steel.forEach(k => materials.push(k));
        this.setMaterials(triCount, materials);
    }

    useGlassyMaterial() {
        var triCount = [];
        for (let i = 0, j = 0; i < this.triangleCount; i += 9, j++) {
            triCount.push(0);
        }
        var materials = [];
        Material.glassy.forEach(k => materials.push(k));
        this.setMaterials(triCount, materials);
    }

    useOneRandomMaterial() {
        var triCount = [];
        for (let i = 0, j = 0; i < this.triangleCount; i += 9, j++) {
            triCount.push(1);
        }
        var materials = [];
        Material.makeRandom().forEach(k => materials.push(k));
        Material.makeRandom().forEach(k => materials.push(k));
        this.setMaterials(triCount, materials);
    }

    setMaterials(triCount, materials) {
        this.textures.materialIndices.image.data.set(triCount);
        this.textures.materials.image.data.set(materials);
        this.textures.materials.needsUpdate = true;
        this.textures.materialIndices.needsUpdate = true;
    }

    createTexture(array, format, type) {
        var texSize = Math.ceil(Math.sqrt(array.length));
        if (texSize < 2048) {
            texSize = Math.pow(2, Math.ceil(Math.log2(texSize)));
        }
        const texArray = new array.constructor(texSize * texSize);
        texArray.set(array);
        const tex = new THREE.DataTexture( texArray, texSize, texSize, format, type );
        tex.flipY = false;
        tex.needsUpdate = true;
        return tex;
    }

    setupRay(p) {
        const { x, y } = p;
        const { width, height } = this.renderer.domElement;

        const uv = new THREE.Vector2((x/width)*2 - 1, -(y/height)*2 + 1);
        uv.x *= width/height;
        
        const direction = new THREE.Vector3( uv.x, uv.y, 1 );
        direction
            .applyMatrix4( this.camera.uvToWorld )
            .sub( this.camera.position )
            .normalize();
        
        const ray = new Ray(this.camera.position.clone(), direction, 0);
        return ray;
    }

    render() {
        if (this.controls.changed || this.frame < parseFloat(window.maxFramesToRender.value)) {

            var dprValue = resAdjust;

            if (this.controls.focusPoint) {
                const ray = this.setupRay({
                    x: this.controls.focusPoint.x * dprValue,
                    y: this.controls.focusPoint.y * dprValue
                });
                this.controls.focusPoint = null;
                const hit = this.voxelGrid.intersect(ray);
                if (hit && hit.obj) {
                    this.focusPoint.copy(add(ray.o, mulS(ray.d, hit.distance)));
                }
            }

            this.focusDistance = this.camera.position.distanceTo(this.focusPoint);

            this.frameTime = Date.now() - this.frameStartTime;
            this.frameStartTime = Date.now();

            const camera = this.camera;
            
            this.material.uniforms.previousCameraMatrixWorld.value.copy(camera.matrixWorld);
            this.material.uniforms.previousCameraProjectionMatrix.value.copy(camera.projectionMatrix);
            this.material.uniforms.previousCameraPosition.value.copy(camera.position);
            camera.setFocalLength(window.focalLength.value);
            camera.lookAt(camera.target);
            camera.updateProjectionMatrix();
            camera.updateMatrixWorld();
            camera.uvToWorld.getInverse(camera.projectionMatrix);
            camera.uvToWorld.multiplyMatrices(camera.matrixWorld, camera.uvToWorld);
            camera.worldToUV.getInverse(camera.matrixWorld);
            camera.worldToUV.multiplyMatrices(camera.worldToUV, camera.projectionMatrix);
            const previousApertureSize = camera.apertureSize;
            camera.apertureSize = 0.01 * window.focalLength.value / window.fStop.value;

            const apertureChanged = (camera.apertureSize !== previousApertureSize)

            const cameraMatrixChanged = !camera.uvToWorld.equals(camera.previousMatrix);

            const controlsActive = (cameraMatrixChanged || apertureChanged);

            if (cameraMatrixChanged || apertureChanged || this.controls.changed || controlsActive) {
                this.frame = 0;
            }

            this.controls.changed = false;

            this.material.uniforms.costVis.value = this.controls.debug;
            this.material.uniforms.aaSize.value = 2;

            this.blitMaterial.uniforms.showConverged.value = !!window.showConverged.checked;
            this.blitMaterial.uniforms.showSampleCount.value = !!window.showSampleCount.checked;

            if (this.renderer.domElement.width !== window.innerWidth*dprValue ||
                this.renderer.domElement.height !== window.innerHeight*dprValue
            ) {
                this.renderer.setSize(window.innerWidth*dprValue, window.innerHeight*dprValue);
                this.renderTarget.setSize(window.innerWidth*dprValue, window.innerHeight*dprValue);
                this.accumRenderTargetA.setSize(window.innerWidth*dprValue, window.innerHeight*dprValue);
                this.accumRenderTargetB.setSize(window.innerWidth*dprValue, window.innerHeight*dprValue);
                this.varianceRenderTarget.setSize(window.innerWidth*dprValue, window.innerHeight*dprValue);
            }

            if (this.frame <= 1) this.rayBudget = this.startRayBudget;

            if (this.frameTime < 18) {
                if (this.frame === 0) this.startRayBudget = Math.min(1000, this.startRayBudget*1.01);
                else this.rayBudget = Math.min(1000, this.rayBudget*1.1);
            } else if (this.frameTime > 20) {
                const clampSpeed = max(0.9, 20 / this.frameTime);
                if (this.frame === 0) this.startRayBudget = Math.max(0.01, this.startRayBudget*clampSpeed);
                else this.rayBudget = Math.max(0.01, this.rayBudget*0.8);
            }

            this.stats.log('Frame', this.frame);
            this.stats.log('Frame time', Math.round(this.frameTime*100)/100 + ' ms');
            this.stats.log('Ray budget', Math.round((this.frame === 0 ? this.startRayBudget : this.rayBudget)*100)/100);

            this.material.uniforms.showFocalPlane.value = !!window.showFocalPlane.checked;
            this.material.uniforms.showBoost.value = !!window.showBoost.checked;
            this.material.uniforms.useTemporalReprojection.value = !!window.useTemporalReprojection.checked;
            this.material.uniforms.temporalReprojectionShowRejection.value = !!window.temporalReprojectionShowRejection.checked;
            this.material.uniforms.temporalReprojectionVarianceCutoff.value = parseFloat(window.temporalReprojectionVarianceCutoff.value);
            this.material.uniforms.temporalReprojectionWeight.value = parseFloat(window.temporalReprojectionWeight.value);
            this.material.uniforms.firstFrameSampleBoost.value = parseFloat(window.firstFrameSampleBoost.value);

            this.material.uniforms.iTime.value = this.frame < 4 ? Math.SQRT2 : (Date.now() - this.startTime) / 1000;
            this.material.uniforms.cameraApertureSize.value = camera.apertureSize;
            this.material.uniforms.iResolution.value[0] = this.renderer.domElement.width;
            this.material.uniforms.iResolution.value[1] = this.renderer.domElement.height;
            this.material.uniforms.rayBudget.value = this.frame === 0 ? this.startRayBudget : this.rayBudget;
            this.material.uniforms.cameraFocusDistance.value = this.focusDistance;

            this.material.uniforms.iFrame.value = this.frame;
            this.blitMaterial.uniforms.iFrame.value = this.frame;
            this.accumMaterial.uniforms.iFrame.value = this.frame;
            this.varianceMaterial.uniforms.iFrame.value = this.frame;

            // Swap render targets for accumulator
            let tmp = this.accumRenderTargetA;
            this.accumRenderTargetA = this.accumRenderTargetB;
            this.accumRenderTargetB = tmp;

            this.material.uniforms.tex.value = this.accumRenderTargetA.texture;
            this.renderer.render(this.scene, camera, this.renderTarget);
            this.accumMaterial.uniforms.accumTex.value = this.accumRenderTargetA.texture;
            this.renderer.render(this.accumMesh, camera, this.accumRenderTargetB);
            this.frame++;

            this.varianceMaterial.uniforms.previousTex.value = this.accumRenderTargetA.texture;
            this.varianceMaterial.uniforms.tex.value = this.accumRenderTargetB.texture;
            this.renderer.render(this.varianceMesh, camera, this.varianceRenderTarget);

            // if (this.frame === 1) {
            //     // Biased start:
            //     //
            //     // 1) Gaussian blur the light
            //     // 2) Use blurred value as starting bias for integration 
            //     //    with weight / kernel width determined by roughness / edginess
            //     //    Rougher => more variance => more weight for blur
            //     //    Edgier => lower blur kernel width

            //     this.material.uniforms.rayBudget.value = 1;

            //     this.material.uniforms.iFrame.value = this.frame;
            //     this.accumMaterial.uniforms.iFrame.value = this.frame;
            //     this.varianceMaterial.uniforms.iFrame.value = this.frame;

            //     // Swap render targets for accumulator
            //     tmp = this.accumRenderTargetA;
            //     this.accumRenderTargetA = this.accumRenderTargetB;
            //     this.accumRenderTargetB = tmp;

            //     this.renderer.render(this.scene, camera, this.renderTarget);
            //     this.accumMaterial.uniforms.accumTex.value = this.accumRenderTargetA.texture;
            //     this.renderer.render(this.accumMesh, camera, this.accumRenderTargetB);
            //     this.frame++;

            //     this.varianceMaterial.uniforms.previousTex.value = this.accumRenderTargetA.texture;
            //     this.varianceMaterial.uniforms.tex.value = this.accumRenderTargetB.texture;
            //     this.renderer.render(this.varianceMesh, camera, this.varianceRenderTarget);

            //     this.material.uniforms.rayBudget.value = 1;

            //     this.material.uniforms.iFrame.value = this.frame;
            //     this.accumMaterial.uniforms.iFrame.value = this.frame;
            //     this.varianceMaterial.uniforms.iFrame.value = this.frame;

            //     // Swap render targets for accumulator
            //     tmp = this.accumRenderTargetA;
            //     this.accumRenderTargetA = this.accumRenderTargetB;
            //     this.accumRenderTargetB = tmp;

            //     this.renderer.render(this.scene, camera, this.renderTarget);
            //     this.accumMaterial.uniforms.accumTex.value = this.accumRenderTargetA.texture;
            //     this.renderer.render(this.accumMesh, camera, this.accumRenderTargetB);
            //     this.frame++;

            //     this.varianceMaterial.uniforms.previousTex.value = this.accumRenderTargetA.texture;
            //     this.varianceMaterial.uniforms.tex.value = this.accumRenderTargetB.texture;
            //     this.renderer.render(this.varianceMesh, camera, this.varianceRenderTarget);
            // }

            if (window.blurMove.checked && this.frame < 50 && !this.controls.debug) {
                this.blurMaterial.uniforms.sigma.value = 5.0 * Math.pow(1.0-(0.5-0.5*Math.cos(Math.PI * this.frame / 50)), 4.0);
                this.blurMaterial.uniforms.direction.value = 0;
                this.blurMaterial.uniforms.tex.value = this.accumRenderTargetB.texture;
                this.renderer.render(this.blurMesh, camera, this.accumRenderTargetA);
                this.blurMaterial.uniforms.direction.value = 1;
                this.blurMaterial.uniforms.tex.value = this.accumRenderTargetA.texture;
                this.renderer.render(this.blurMesh, camera, this.renderTarget);
                this.blitMaterial.uniforms.tex.value = this.renderTarget.texture;
            } else {
                this.bloomMaterial.uniforms.tex.value = this.accumRenderTargetB.texture;
                this.renderer.render(this.bloomMesh, camera, this.renderTarget);
                this.blitMaterial.uniforms.tex.value = this.accumRenderTargetB.texture;
            }
            this.renderer.render(this.blitMesh, camera);

            camera.previousMatrix.copy(camera.uvToWorld);

            if (this.frame === 0 && !this.firstFrameTimerFired) {
                console.timeEnd('Load to first frame');
                this.firstFrameTimerFired = true;
            }
        }
    }
}


function LoadOBJ(path) {
    return new Promise((resolve, reject) => {
        new THREE.OBJLoader().load(path, resolve, null, reject);
    });
};

(async function() {
    console.time('Load to first frame');

    var tracer = false;
    var voxelGrid = false;
    var shaders = false;

    var loadCount = 0;
    const onLoad = function() {
        if (loadCount === 4) {
            if (!voxelGrid || !shaders) {
                setTimeout(onLoad, 10);
                return;
            }
            tracer = new WebGLTracer2(voxelGrid, shaders.join('\n'), blueNoiseTexture, diffuseTexture, metallicTexture, roughnessTexture, normalTexture);

            const tick = () => {
                tracer.render();
                requestAnimationFrame(tick);
            }
            tick();
        
            document.body.append(tracer.renderer.domElement);
        } else {
            loadCount++;
        }
    };

    const blueNoiseTexture = new THREE.TextureLoader().load('blue_noise.png', onLoad);
    blueNoiseTexture.wrapS = THREE.RepeatWrapping;
    blueNoiseTexture.wrapT = THREE.RepeatWrapping;
    const diffuseTexture = new THREE.TextureLoader().load('Metal_basecolor.png', onLoad);
    diffuseTexture.wrapS = THREE.RepeatWrapping;
    diffuseTexture.wrapT = THREE.RepeatWrapping;
    const metallicTexture = new THREE.TextureLoader().load('Metal_metallic.png', onLoad);
    metallicTexture.wrapS = THREE.RepeatWrapping;
    metallicTexture.wrapT = THREE.RepeatWrapping;
    const normalTexture = new THREE.TextureLoader().load('Metal_normal.png', onLoad);
    normalTexture.wrapS = THREE.RepeatWrapping;
    normalTexture.wrapT = THREE.RepeatWrapping;
    const roughnessTexture = new THREE.TextureLoader().load('Metal_roughness.png', onLoad);
    roughnessTexture.wrapS = THREE.RepeatWrapping;
    roughnessTexture.wrapT = THREE.RepeatWrapping;

    const shaderNames = ['primitives', 'voxelgrid_superflat', 'trace2'];

    const shaderRes = shaderNames.map(name => fetch(`lib/${name}.glsl`));
    const bunny = await ObjParse.load('bunny.obj');

    const shadersT = await Promise.all(shaderRes.map(async res => (await res).text()));
    shaders = shadersT;

    console.time('OBJ munging');
    var verts = bunny.vertices;
    var normals = bunny.normals;

    var bbox = { 
        min: vec3(Infinity),
        max: vec3(-Infinity)
    };
    for (let i = 0; i < verts.length; i += 3) {
        const x = verts[i];
        const y = verts[i+1];
        const z = verts[i+2];
        if (x < bbox.min.x) bbox.min.x = x;
        if (x > bbox.max.x) bbox.max.x = x;
        if (y < bbox.min.y) bbox.min.y = y;
        if (y > bbox.max.y) bbox.max.y = y;
        if (z < bbox.min.z) bbox.min.z = z;
        if (z > bbox.max.z) bbox.max.z = z;
    }

    var scale = 2.5 / max(bbox.max.z - bbox.min.z, bbox.max.x - bbox.min.x);
    var xOffset = -(bbox.max.x+bbox.min.x)/2;
    var zOffset = -(bbox.max.z+bbox.min.z)/2;
    var yOffset = -bbox.min.y;

    for (var i = 0; i < verts.length; i += 3) {
        verts[i] += xOffset;
        verts[i] *= scale;
        verts[i+1] += yOffset;
        verts[i+1] *= scale;
        verts[i+2] += zOffset;
        verts[i+2] *= scale;
    }
    bbox.min = mulS(add(bbox.min, vec3(xOffset, yOffset, zOffset)), scale);
    bbox.max = mulS(add(bbox.max, vec3(xOffset, yOffset, zOffset)), scale);

    var bunnyTris = [];

    var color = vec3(0.85, 0.53, 0.15);
    for (var i = 0; i < verts.length; i += 3*3) {
        var u = vec3(verts[i], verts[i+1], verts[i+2]);
        var v = vec3(verts[i+3], verts[i+4], verts[i+5]);
        var w = vec3(verts[i+6], verts[i+7], verts[i+8]);
        var x = vec3(normals[i], normals[i+1], normals[i+2]);
        var y = vec3(normals[i+3], normals[i+4], normals[i+5]);
        var z = vec3(normals[i+6], normals[i+7], normals[i+8]);
        bunnyTris.push(new Triangle([u,v,w], [x,y,z], color));
    }
    bunnyTris.bbox = bbox;

    var size = sub(bunnyTris.bbox.max, bunnyTris.bbox.min);
    var m = Math.max(size.x, size.y, size.z);
    // decent for dragon [32, 4, 2, 2] and [16,4,4,2] and [32,2,2,2]
    var grid = [16,16]; //[32,4,2,2];
    if (bunnyTris.length < 10000) {
        // Use low-res grid
        // Fastest JS exec: [32]
        // Nice mix of VG steps + intersects: [4,4,4]
        // + Fast JS exec: [8, 8]
        grid = [8,8];
    }
    console.timeEnd('OBJ munging');

    console.time('Create VoxelGrid');
    const vg  = new VoxelGrid4(bunnyTris.bbox.min, vec3(m), grid, 0);
    vg.addTriangles(bunnyTris);
    voxelGrid = vg;
    console.timeEnd('Create VoxelGrid');

    document.querySelectorAll("input[type=range]").forEach(el => {
        el.oninput = function(ev) {
            tracer.controls.pinching = true;
            tracer.controls.changed = true;
            this.setAttribute('value', this.value);
        };
        el.onchange = function(ev) {
            tracer.controls.pinching = false;
            tracer.controls.changed = true;
            this.setAttribute('value', this.value);
        };
    });
    

    window.focalLength.oninput = function(ev) {
        tracer.controls.pinching = true;
        tracer.controls.changed = true;
        var f = this.value / tracer.camera.getFocalLength();
        tracer.controls.distance *= f;
        tracer.controls.updateCameraPosition();
        this.setAttribute('value', this.value);
    };    
    window.focalLength.onchange = function(ev) {
        tracer.controls.pinching = false;
        tracer.controls.changed = true;
        var f = this.value / tracer.camera.getFocalLength();
        tracer.controls.distance *= f;
        tracer.controls.updateCameraPosition();
        this.setAttribute('value', this.value);
    };
        
    window.showFocalPlane.onchange = 
    window.showBoost.onchange = 
    window.blurMove.onchange = function() {
        tracer.controls.changed = true;
    };

    window.useRandomMaterials.onclick = function() {
        tracer.useRandomMaterials();
        tracer.controls.changed = true;
    };

    window.useSteelMaterial.onclick = function() {
        tracer.useSteelMaterial();
        tracer.controls.changed = true;
    };

    window.useGlassyMaterial.onclick = function() {
        tracer.useGlassyMaterial();
        tracer.controls.changed = true;
    };

    window.useOneRandomMaterial.onclick = function() {
        tracer.useOneRandomMaterial();
        tracer.controls.changed = true;
    };

})();
