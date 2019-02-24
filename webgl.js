const dpr = window.devicePixelRatio || 1;

class WebGLTracer {
    constructor(vgArray, traceGLSL) {
        var canvas = document.createElement( 'canvas' );
        var context = canvas.getContext( 'webgl2' );
        const texSize = 1024;
        const paddedVgArray = new Float32Array(texSize * texSize);
        paddedVgArray.set(vgArray);
        this.vgArray = paddedVgArray;
        this.vgTexture = new THREE.DataTexture( paddedVgArray, texSize, texSize, THREE.RedFormat, THREE.FloatType );
        this.vgTexture.flipY = false;
        this.vgTexture.needsUpdate = true;
        this.ivgTexture = new THREE.DataTexture( new Int32Array( paddedVgArray.buffer ), texSize, texSize, THREE.RedIntegerFormat, THREE.IntType );
        this.ivgTexture.flipY = false;
        this.ivgTexture.needsUpdate = true;
        this.renderer = new THREE.WebGLRenderer({ canvas, context });
        this.renderer.setSize(window.innerWidth*dpr, window.innerHeight*dpr);
        this.renderer.setClearColor(0x00ffff);
        this.renderer.clear();

        var camera = new THREE.PerspectiveCamera(55, 1, 0.1, 100);
        camera.target = new THREE.Vector3(0, 2.75, 0);
        camera.focusPoint = vec3(-1, 0.7, 0.2);
		camera.apertureSize = Math.pow(1.33, -3);
        camera.inverseMatrix = new THREE.Matrix4()
        camera.lookAt(camera.target);
        camera.updateProjectionMatrix();
        camera.updateMatrixWorld();

        this.material = new THREE.ShaderMaterial({
            uniforms: {
                iTime: { value: 1.0 },
                arrayTex: { value: this.vgTexture },
                iarrayTex: { value: this.ivgTexture },
                arrayTexWidth: { value: texSize },
                iResolution: { value: [canvas.width, canvas.height] },
                cameraFocusPoint: { value: camera.focusPoint },
                focusDistance: { value: 1.0 },
                cameraApertureSize: { value: 0.3 },
                cameraInverseMatrix: { value: camera.inverseMatrix }
            },
            vertexShader: `
            #version 300 es

            precision highp float;
            precision highp int;

            void main() {
                gl_Position = vec4(position.xyz, 1.0);
            }
            `,
            fragmentShader: traceGLSL + `

            out vec4 FragColor;    

            #define AA_SIZE 1.0

            void main() {
                Array array = Array(arrayTexWidth);
                vec3 sum = vec3(0.0);
                for (float y = 0.0; y < AA_SIZE; y++)
                for (float x = 0.0; x < AA_SIZE; x++) {
                    vec3 c = trace(array, gl_FragCoord.xy + vec2(x,y) / AA_SIZE);
                    c = -exp(-c * 0.5) + 1.0;
                    sum +=  c;
                }
                FragColor = vec4(sum / (AA_SIZE*AA_SIZE), 1.0);
            }
            `
        });
        this.mesh = new THREE.Mesh(
            new THREE.PlaneGeometry(2,2,1),
            this.material
        );
        this.material.depthTest = false;
        this.material.depthWrite = false;
        this.mesh.frustumCulled = false;
        this.mesh.position.z = -0.5;
        this.mesh.rotation.y = Math.PI;

        this.scene = new THREE.Scene();
        this.camera = camera;

        this.scene.add(this.camera);
        this.scene.add(this.mesh);

        this.controls = new CameraControls(camera, canvas);
        camera.positionOffset.y = 0.1;
        this.startTime = Date.now();

        window.onresize = () => {
            this.renderer.setSize(window.innerWidth*dpr, window.innerHeight*dpr);
        }
    }

    render() {
        const camera = this.camera;
        camera.lookAt(camera.target);
		camera.updateProjectionMatrix();
        camera.updateMatrixWorld();
        camera.inverseMatrix.getInverse(camera.projectionMatrix);
        camera.inverseMatrix.multiplyMatrices(camera.matrixWorld, camera.inverseMatrix);
        this.material.uniforms.iTime.value = (Date.now() - this.startTime) / 1000;
        this.material.uniforms.iResolution.value[0] = this.renderer.domElement.width;
        this.material.uniforms.iResolution.value[1] = this.renderer.domElement.height;

        this.renderer.render(this.scene, camera);
    }
}

(async function() {
    const bunnyRes = await fetch('bunny.vg3');
    const vgRes = await fetch('lib/voxelgrid.glsl');
    const traceRes = await fetch('lib/trace.glsl');
    const vgText = await vgRes.text();
    const traceText = await traceRes.text();
    const bunnyVG = await bunnyRes.arrayBuffer();

    const tracer = new WebGLTracer(new Float32Array(bunnyVG), vgText + '\n' + traceText);
    tracer.render();

    const tick = () => {

        // tracer.camera.position.x = Math.cos(Date.now()/1000) * 15;
        // tracer.camera.position.y = 9.0 + Math.cos(Date.now()/4531) * 2.4;
        // tracer.camera.position.z = Math.sin(Date.now()/1000) * 15;
        // tracer.camera.fov = 120;
        tracer.render();
        requestAnimationFrame(tick);
    }
    tick();

    document.body.append(tracer.renderer.domElement);
})();
