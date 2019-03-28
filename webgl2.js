const mobile = /mobile/i.test(navigator.userAgent);
const dpr = (window.devicePixelRatio || 1);

class WebGLTracer2 {
    constructor(vg, traceGLSL, blueNoiseTexture) {
        var canvas = document.createElement( 'canvas' );
        var context = canvas.getContext( 'webgl2' );

        this.voxelGrid = vg;

        console.time('Serialize VoxelGrid');
        const arrays = vg.serialize();
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

        this.textures = {
            voxelIndex: this.createTexture(new Int16Array(arrays.voxelIndex.buffer), THREE.RedIntegerFormat, THREE.ShortType),
            triIndices: this.createTexture(arrays.triIndices, THREE.RedIntegerFormat, THREE.UnsignedShortType),
            triangles: this.createTexture(arrays.triangles, THREE.RedFormat, THREE.FloatType),
            normals: this.createTexture(arrays.normals, THREE.RedFormat, THREE.FloatType)
        };

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
        camera.inverseMatrix = new THREE.Matrix4()
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

                varianceTexture: { value: this.varianceRenderTarget.texture },
                
                vgOrigin: { value: new THREE.Vector3().copy(vg._origin) },
                vgScale: { value: vg._scale },
                vgSize: { value: vg._size },

                blueNoise: { value: blueNoiseTexture },

                iResolution: { value: [canvas.width, canvas.height] },
                cameraFocusDistance: { value: this.focusDistance },
                cameraApertureSize: { value: camera.apertureSize },
                cameraInverseMatrix: { value: camera.inverseMatrix },
                cameraMatrixWorld: { value: camera.matrixWorld },
                deviceEpsilon: {value: mobile ? 0.01 : 0.0001},
                deviceEpsilonTrace: {value: mobile ? 0.05 : 0.01},
                roughness: {value: 0.2},
                costVis: {value: false},
                aaSize: {value: 4.0},
                stripes: { value: false },
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
            uniform sampler2D varianceTexture;
            
            ${traceGLSL}

            uniform float aaSize;
            uniform float rayBudget;
            uniform bool showBoost;

            out vec4 FragColor;    

            void main() {

                float samples = 1.0;
                vec4 sum = vec4(0.0);

                if (iFrame > 1.0) {

                    vec4 varianceMetrics = texelFetch(varianceTexture, ivec2(gl_FragCoord.xy), 0);
                    float errorLum = varianceMetrics.x;
                    float totalVariance = varianceMetrics.y;
                    bool converged = varianceMetrics.z > 0.0;
                    bool convergedVariance = varianceMetrics.w > 0.0;

                    if (rayBudget < 4.0 && ((iFrame > 2.0 && iFrame < 10.0) || (iFrame > 100.0 && errorLum < 0.0001)) && (converged || convergedVariance || (iFrame >= 2.0 && (errorLum < 0.01 || totalVariance < 0.01)))) {
                        FragColor = vec4(0.0, 0.0, 0.0, 0.0);
                        return;
                    }

                    vec2 centerToPixel = vec2(iResolution.x / iResolution.y, 1.0) * (gl_FragCoord.xy / iResolution.xy - 0.5);
                    vec2 centerToPixel8x8 = vec2(iResolution.x / iResolution.y, 1.0) * (8.0*floor(gl_FragCoord.xy / 8.0) / iResolution.xy - 0.5);
                    float distanceToCenterPx = length(centerToPixel);
                    float distanceToCenter = length(centerToPixel8x8);
                    float boostRadius = max(0.2, rayBudget * 0.4);
                    float fullBoostRadius = max(0.025, rayBudget * 0.1);
                    float boostExponent = 8.0 / max(0.125, rayBudget);
                    float fullBoostSamples = rayBudget * 5.0;
                    float boostFactor = (1.0 - clamp((distanceToCenter - fullBoostRadius) / boostRadius, 0.0, 1.0));
                    boostFactor += totalVariance/4.0;
                    boostFactor = clamp(boostFactor, 0.0, 1.0);
                    float boost = pow(boostFactor, boostExponent) * fullBoostSamples;
                    float sampleCountJitter = 0.0; //fract(random(vec2(((1.0+sqrt(2.0)))+iTime, (9.0+sqrt(221.0))*0.1+iTime)));

                    samples = boost + sampleCountJitter;

                    if (rayBudget > 4.0 && iFrame > 100.0) {
                        samples = 1.0 + totalVariance;
                    }
                    
                    if (showBoost) {
                        sum.g += samples;
                        if (distanceToCenterPx < fullBoostRadius) {
                            if (distanceToCenterPx < 0.002 || (fullBoostRadius - distanceToCenterPx) < 0.002) {
                                sum.r += samples;
                                sum.gb *= 0.0;
                            }
                        }
                    }
                }

                for (float i = 0.0; i < samples; i++) {
                    float y = mod(i, 2.0);
                    float x = i - y * 2.0;
                    float ry = mod(y + iFrame / aaSize, aaSize);
                    float rx = mod(x + iFrame - aaSize * ry, aaSize);
                    vec3 c = trace(gl_FragCoord.xy + (vec2(rx,ry) + vec2(random(i+0.5*vec2(rx,ry)), random(i+0.5+0.5*vec2(rx,ry)))) / aaSize);
                    sum += vec4(c, 1.0);
                }
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

                if (iFrame == 0.0) {
                    FragColor = vec4(1.0, 1.0, 0.0, 0.0);
                    return;
                }

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

        this.blitMaterial = new THREE.RawShaderMaterial({
            uniforms: {
                tex: { value: this.renderTarget.texture },
                varianceTexture: { value: this.varianceRenderTarget.texture },
                showConverged: { value: false },
                showSampleCount: { value: false }
            },
            vertexShader: this.material.vertexShader,
            fragmentShader: `#version 300 es

            precision highp float;
            precision highp int;
            
            uniform sampler2D tex;
            uniform sampler2D varianceTexture;

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

                if (errorLum > 0.01 || length(FragColor.rgb) < 0.4) {
                    vec3 v[9];
                    v[0] = toGamma(texelFetch(tex, ivec2(gl_FragCoord.xy) + ivec2(0, -1), 0));
                    v[1] = toGamma(texelFetch(tex, ivec2(gl_FragCoord.xy) + ivec2(0, -1), 0));
                    v[2] = toGamma(texelFetch(tex, ivec2(gl_FragCoord.xy) + ivec2(1, 0), 0));
                    v[3] = toGamma(texelFetch(tex, ivec2(gl_FragCoord.xy) + ivec2(-1, 0), 0));
                    v[4] = FragColor.rgb;
                    v[5] = toGamma(texelFetch(tex, ivec2(gl_FragCoord.xy) + ivec2(1, 0), 0));
                    v[6] = toGamma(texelFetch(tex, ivec2(gl_FragCoord.xy) + ivec2(0, 1), 0));
                    v[7] = toGamma(texelFetch(tex, ivec2(gl_FragCoord.xy) + ivec2(0, 1), 0));
                    v[8] = toGamma(texelFetch(tex, ivec2(gl_FragCoord.xy) + ivec2(1, 0), 0));

                    if (length(v[4]) > 0.2 && (length(v[4]) > 1.2 || length(v[4] - v[5]) > 0.1 || length(v[4] - v[3]) > 0.1 || length(v[4] - v[1]) > 0.1 || length(v[4] - v[7]) > 0.1)) { 
                        #define s2(a, b)				temp = a; a = min(a, b); b = max(temp, b);
                        #define mn3(a, b, c)			s2(a, b); s2(a, c);
                        #define mx3(a, b, c)			s2(b, c); s2(a, c);

                        #define mnmx3(a, b, c)			mx3(a, b, c); s2(a, b);                                   // 3 exchanges
                        #define mnmx4(a, b, c, d)		s2(a, b); s2(c, d); s2(a, c); s2(b, d);                   // 4 exchanges
                        #define mnmx5(a, b, c, d, e)	s2(a, b); s2(c, d); mn3(a, c, e); mx3(b, d, e);           // 6 exchanges
                        #define mnmx6(a, b, c, d, e, f) s2(a, d); s2(b, e); s2(c, f); mn3(a, b, c); mx3(d, e, f); // 7 exchanges

                        vec3 temp;
                        mnmx6(v[0], v[1], v[2], v[3], v[4], v[5]);
                        mnmx5(v[1], v[2], v[3], v[4], v[6]);
                        mnmx4(v[2], v[3], v[4], v[7]);
                        mnmx3(v[3], v[4], v[8]);
                        FragColor.rgb = v[4];
                    } else {
                        FragColor.rgb = toGamma(
                            src +
                            texelFetch(tex, ivec2(gl_FragCoord.xy) + ivec2(0, -1), 0) + 
                            texelFetch(tex, ivec2(gl_FragCoord.xy) + ivec2(0, 1), 0) + 
                            texelFetch(tex, ivec2(gl_FragCoord.xy) + ivec2(1, 0), 0) + 
                            texelFetch(tex, ivec2(gl_FragCoord.xy) + ivec2(-1, 0), 0)
                        );
                    }
                }

                if (showConverged) {
                    if ((converged || convergedVariance || errorLum < 0.01 || totalVariance < 0.01 )) {
                        FragColor.r -= 0.5;
                        FragColor.b += 0.5;
                    } else {
                        FragColor.r += 0.5;
                        FragColor.b -= 0.5;
                    }
                }
                if (showSampleCount) {
                    FragColor = vec4(src.aaa / 100.0, 1.0);
                }
                // FragColor = vec4(errorLum/2.0, 0.0, src.a/10.0, 1.0);
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
        this.movingRayBudget = 1;
        this.frameTime = 0;
        this.frameStartTime = Date.now();

        window.onresize = () => {
            this.controls.changed = true;
        };
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
        
        const direction = new THREE.Vector3( uv.x, uv.y, 1 );
        direction
            .applyMatrix4( this.camera.inverseMatrix )
            .sub( this.camera.position )
            .normalize();
        
        const ray = new Ray(this.camera.position.clone(), direction, 0);
        return ray;
    }

    render() {
        if (true) {

            if (this.controls.focusPoint) {
                const ray = this.setupRay(this.controls.focusPoint);
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
            camera.lookAt(camera.target);
            camera.updateProjectionMatrix();
            camera.updateMatrixWorld();
            camera.inverseMatrix.getInverse(camera.projectionMatrix);
            camera.inverseMatrix.multiplyMatrices(camera.matrixWorld, camera.inverseMatrix);
            const previousApertureSize = camera.apertureSize;
            camera.apertureSize = Math.pow(1.33, 1-window.apertureSize.value);

            const apertureChanged = (camera.apertureSize !== previousApertureSize)

            const cameraMatrixChanged = !camera.inverseMatrix.equals(camera.previousMatrix);

            const controlsActive = (cameraMatrixChanged || apertureChanged || this.controls.down || this.controls.pinching);

            if (cameraMatrixChanged || apertureChanged || this.controls.changed) {
                this.rayBudget = 1;
                this.frame = 0;
            }

            if (this.frameTime < 20) {
                this.rayBudget = Math.min(1000, this.rayBudget*1.05);
            } else if (this.frameTime > 40) {
                this.rayBudget = Math.max(0.01, this.rayBudget*0.8);
            }
            this.stats.log('Frame time', Math.round(this.frameTime*100)/100 + ' ms');
            this.stats.log('Ray budget', Math.round(this.rayBudget*100)/100);
            
            this.controls.changed = false;

            this.material.uniforms.costVis.value = this.controls.debug;
            this.material.uniforms.aaSize.value = 2;

            this.blitMaterial.uniforms.showConverged.value = !!window.showConverged.checked;
            this.blitMaterial.uniforms.showSampleCount.value = !!window.showSampleCount.checked;

            var dprValue = dpr;

            if (controlsActive) {
                this.rayBudget = 1;
                if (this.renderer.domElement.width !== window.innerWidth ||
                    this.renderer.domElement.height !== window.innerHeight
                ) {
                    this.rayBudget = 1;
                    this.frame = 0;
                    this.renderer.setSize(window.innerWidth, window.innerHeight);
                    this.renderTarget.setSize(window.innerWidth, window.innerHeight);
                    this.accumRenderTargetA.setSize(window.innerWidth, window.innerHeight);
                    this.accumRenderTargetB.setSize(window.innerWidth, window.innerHeight);
                    this.varianceRenderTarget.setSize(window.innerWidth, window.innerHeight);
                }
                dprValue = 1;
            } else {
                if (this.renderer.domElement.width !== window.innerWidth*dpr ||
                    this.renderer.domElement.height !== window.innerHeight*dpr
                ) {
                    this.rayBudget /= dpr * dpr;
                    this.frame = 0;
                    this.renderer.setSize(window.innerWidth*dpr, window.innerHeight*dpr);
                    this.renderTarget.setSize(window.innerWidth*dpr, window.innerHeight*dpr);
                    this.accumRenderTargetA.setSize(window.innerWidth*dpr, window.innerHeight*dpr);
                    this.accumRenderTargetB.setSize(window.innerWidth*dpr, window.innerHeight*dpr);
                    this.varianceRenderTarget.setSize(window.innerWidth*dpr, window.innerHeight*dpr);
                }
            }

            this.material.uniforms.stripes.value = !!window.stripes.checked;
            this.material.uniforms.showFocalPlane.value = !!window.showFocalPlane.checked;
            this.material.uniforms.showBoost.value = !!window.showBoost.checked;

            this.material.uniforms.iTime.value = this.frame < 4 ? Math.SQRT2 : (Date.now() - this.startTime) / 1000;
            this.material.uniforms.cameraApertureSize.value = camera.apertureSize;
            this.material.uniforms.iResolution.value[0] = this.renderer.domElement.width;
            this.material.uniforms.iResolution.value[1] = this.renderer.domElement.height;
            this.material.uniforms.roughness.value = window.roughness.value / 100;
            this.material.uniforms.rayBudget.value = this.rayBudget;
            this.material.uniforms.cameraFocusDistance.value = this.focusDistance;

            this.material.uniforms.iFrame.value = this.frame;
            this.accumMaterial.uniforms.iFrame.value = this.frame;
            this.varianceMaterial.uniforms.iFrame.value = this.frame;

            // Swap render targets for accumulator
            let tmp = this.accumRenderTargetA;
            this.accumRenderTargetA = this.accumRenderTargetB;
            this.accumRenderTargetB = tmp;

            this.renderer.render(this.scene, camera, this.renderTarget);
            this.accumMaterial.uniforms.accumTex.value = this.accumRenderTargetA.texture;
            this.renderer.render(this.accumMesh, camera, this.accumRenderTargetB);
            this.frame++;

            this.varianceMaterial.uniforms.previousTex.value = this.accumRenderTargetA.texture;
            this.varianceMaterial.uniforms.tex.value = this.accumRenderTargetB.texture;
            this.renderer.render(this.varianceMesh, camera, this.varianceRenderTarget);

            if (this.frame === 1) {
                this.material.uniforms.rayBudget.value = 1;

                this.material.uniforms.iFrame.value = this.frame;
                this.accumMaterial.uniforms.iFrame.value = this.frame;
                this.varianceMaterial.uniforms.iFrame.value = this.frame;

                // Swap render targets for accumulator
                tmp = this.accumRenderTargetA;
                this.accumRenderTargetA = this.accumRenderTargetB;
                this.accumRenderTargetB = tmp;

                this.renderer.render(this.scene, camera, this.renderTarget);
                this.accumMaterial.uniforms.accumTex.value = this.accumRenderTargetA.texture;
                this.renderer.render(this.accumMesh, camera, this.accumRenderTargetB);
                this.frame++;

                this.varianceMaterial.uniforms.previousTex.value = this.accumRenderTargetA.texture;
                this.varianceMaterial.uniforms.tex.value = this.accumRenderTargetB.texture;
                this.renderer.render(this.varianceMesh, camera, this.varianceRenderTarget);

                this.material.uniforms.rayBudget.value = 1;

                this.material.uniforms.iFrame.value = this.frame;
                this.accumMaterial.uniforms.iFrame.value = this.frame;
                this.varianceMaterial.uniforms.iFrame.value = this.frame;

                // Swap render targets for accumulator
                tmp = this.accumRenderTargetA;
                this.accumRenderTargetA = this.accumRenderTargetB;
                this.accumRenderTargetB = tmp;

                this.renderer.render(this.scene, camera, this.renderTarget);
                this.accumMaterial.uniforms.accumTex.value = this.accumRenderTargetA.texture;
                this.renderer.render(this.accumMesh, camera, this.accumRenderTargetB);
                this.frame++;

                this.varianceMaterial.uniforms.previousTex.value = this.accumRenderTargetA.texture;
                this.varianceMaterial.uniforms.tex.value = this.accumRenderTargetB.texture;
                this.renderer.render(this.varianceMesh, camera, this.varianceRenderTarget);
            }

            if (window.blurMove.checked && this.frame < 30 && dprValue === 1 && !this.controls.debug) {
                this.blurMaterial.uniforms.sigma.value = 15.0 * Math.pow(1.0-(0.5-0.5*Math.cos(Math.PI * this.frame / 30)), 4.0);
                this.blurMaterial.uniforms.direction.value = 0;
                this.blurMaterial.uniforms.tex.value = this.accumRenderTargetB.texture;
                this.renderer.render(this.blurMesh, camera, this.accumRenderTargetA);
                this.blurMaterial.uniforms.direction.value = 1;
                this.blurMaterial.uniforms.tex.value = this.accumRenderTargetA.texture;
                this.renderer.render(this.blurMesh, camera, this.renderTarget);
                this.blitMaterial.uniforms.tex.value = this.renderTarget.texture;
            } else {
                this.blitMaterial.uniforms.tex.value = this.accumRenderTargetB.texture;
            }
            this.renderer.render(this.blitMesh, camera);

            camera.previousMatrix.copy(camera.inverseMatrix);

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

    const onLoad = function() {
        if (tracer) {
            tracer.render();

            const tick = () => {
                tracer.render();
                requestAnimationFrame(tick);
            }
            tick();
        
            document.body.append(tracer.renderer.domElement);
        } else {
            setTimeout(onLoad, 10);
        }
    };

    const blueNoiseTexture = new THREE.TextureLoader().load('blue_noise.png', onLoad);

    const shaderNames = ['primitives', 'voxelgrid_superflat', 'trace2'];

    const shaderRes = shaderNames.map(name => fetch(`lib/${name}.glsl`));
    const bunny = await ObjParse.load('bunny.obj');

    const shaders = await Promise.all(shaderRes.map(async res => (await res).text()));

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

    var scale = 2.0 / (bbox.max.y - bbox.min.y);
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
        // grid = [8,8];
    }
    console.timeEnd('OBJ munging');

    console.time('Create VoxelGrid');
    const voxelGrid = new VoxelGrid4(bunnyTris.bbox.min, vec3(m), grid, 0);
    voxelGrid.addTriangles(bunnyTris);
    console.timeEnd('Create VoxelGrid');


    tracer = new WebGLTracer2(voxelGrid, shaders.join('\n'), blueNoiseTexture);

    
    window.roughness.oninput = window.apertureSize.oninput = function(ev) {
        tracer.controls.pinching = true;
        tracer.controls.changed = true;
        this.setAttribute('value', this.value);
    };

    window.roughness.onchange = window.apertureSize.onchange = function() {
        tracer.controls.pinching = false;
        tracer.controls.changed = true;
        this.setAttribute('value', this.value);
    };
        
    window.showFocalPlane.onchange = 
    window.showBoost.onchange = 
    window.stripes.onchange = 
    window.blurMove.onchange = function() {
        tracer.controls.changed = true;
    };

})();
