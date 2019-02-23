const dpr = window.devicePixelRatio || 1;

class WebGLTracer {
    constructor(vgArray, traceGLSL) {
        var canvas = document.createElement( 'canvas' );
        var context = canvas.getContext( 'webgl2' );
        const paddedVgArray = new Float32Array(8192 * 8192);
        paddedVgArray.set(vgArray);
        this.vgArray = paddedVgArray;
        // for (var i=0; i<900; i++) { 
        //     paddedVgArray[i*3 + 0] = 1.0;
        //     paddedVgArray[i*3 + 1] = 0.0;
        //     paddedVgArray[i*3 + 2] = 1.0;    
        // }
        this.vgTexture = new THREE.DataTexture( paddedVgArray, 1024, 1024, THREE.RedFormat, THREE.FloatType );
        this.vgTexture.flipY = false;
        this.vgTexture.needsUpdate = true;
        this.renderer = new THREE.WebGLRenderer({ canvas, context });
        this.renderer.setSize(window.innerWidth*dpr, window.innerHeight*dpr);
        this.renderer.setClearColor(0x00ffff);
        this.renderer.clear();

        var camera = new THREE.PerspectiveCamera(55, 1, 0.1, 100);
        camera.target = new THREE.Vector3(0, 0.75, 0);
        camera.focusPoint = vec3(-1, 0.7, 0.2);
		camera.apertureSize = Math.pow(1.33, -3);
        camera.inverseProjectionMatrix = new THREE.Matrix4()
        camera.lookAt(camera.target);
        camera.updateProjectionMatrix();
        camera.updateMatrixWorld();
        camera.inverseProjectionMatrix.getInverse(camera.projectionMatrix);        

        this.material = new THREE.ShaderMaterial({
            uniforms: {
                iTime: { value: 1.0 },
                arrayTex: { value: this.vgTexture },
                arrayTexWidth: { value: 1024 },
                iResolution: { value: [window.innerWidth, window.innerHeight] },
                cameraFocusPoint: { value: camera.focusPoint },
                focusDistance: { value: 1.0 },
                cameraApertureSize: { value: 0.3 },
                cameraMatrixWorld: { value: camera.matrixWorld },
                cameraInverseProjectionMatrix: { value: camera.inverseProjectionMatrix }
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

            void main() {
                Array array = Array(arrayTexWidth);
                vec3 c = trace(array);
                c = -exp(-c * 0.5) + 1.0;
                FragColor = vec4(c, 1.0);
            }
            `
        });
        this.mesh = new THREE.Mesh(
            new THREE.PlaneGeometry(2,2,1),
            this.material
        );
        this.mesh.position.z = -0.5;
        this.mesh.rotation.y = Math.PI;

        this.scene = new THREE.Scene();
        this.camera = camera;

        this.scene.add(this.camera);
        this.scene.add(this.mesh);

        this.controls = new CameraControls(camera, canvas);
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
        camera.inverseProjectionMatrix.getInverse(camera.projectionMatrix);
        this.material.uniforms.iTime.value = (Date.now() - this.startTime) / 1000;
        this.material.uniforms.iResolution.value[0] = window.innerWidth*dpr;
        this.material.uniforms.iResolution.value[1] = window.innerHeight*dpr;

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
        tracer.render();
        requestAnimationFrame(tick);
    }
    tick();

    document.body.append(tracer.renderer.domElement);
})();
