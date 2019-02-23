class WebGLTracer {
    constructor(vgArray, traceGLSL) {
        var canvas = document.createElement( 'canvas' );
        var context = canvas.getContext( 'webgl2' );
        this.vgArray = vgArray;
        this.vgTexture = new THREE.DataTexture( vgArray, 2048, Math.ceil(vgArray.length / 2048), THREE.LuminosityType, THREE.FloatType );
        this.renderer = new THREE.WebGLRenderer({ canvas, context });
        this.renderer.setSize(300, 300);
        this.renderer.setClearColor(0x00ffff);
        this.renderer.clear();

        this.material = new THREE.ShaderMaterial({
            uniforms: {
                time: { value: 1.0 },
                arrayTex: { value: this.vgTexture },
                arrayTexWidth: { value: 2048 },
                iResolution: { value: new THREE.Vector2(this.renderer.width, this.renderer.height) },
                cameraFocusPoint: { value: new THREE.Vector3() },
                focusDistance: { value: 1.0 },
                cameraApertureSize: { value: 0.3 }
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

        var camera = new THREE.PerspectiveCamera(55, 1, 0.1, 100);
        camera.target = new THREE.Vector3(0, 0.75, 0);
        camera.focusPoint = vec3(-1, 0.7, 0.2);
        camera.positionOffset = new THREE.Vector3();
		camera.apertureSize = Math.pow(1.33, -3);
        camera.lookAt(camera.target);
        camera.updateProjectionMatrix();
        camera.updateMatrixWorld();
    
        this.scene = new THREE.Scene();
        this.camera = camera;

        this.scene.add(this.camera);
        this.scene.add(this.mesh);
    }

    render() {
        const camera = this.camera;
        camera.lookAt(camera.target);
		camera.updateProjectionMatrix();
        camera.updateMatrixWorld();

        this.renderer.render(this.scene, this.camera);
    }
}

(async function() {
    const bunnyRes = await fetch('bunny.vg3');
    const vgRes = await fetch('lib/voxelgrid.glsl');
    const traceRes = await fetch('lib/trace.glsl');
    const vgText = await vgRes.text();
    const traceText = await traceRes.text();
    const bunnyVG = await bunnyRes.arrayBuffer();

    const tracer = new WebGLTracer(bunnyVG, vgText + '\n' + traceText);
    tracer.render();

    document.body.append(tracer.renderer.domElement);
})();
