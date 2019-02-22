class WebGLTracer {
    constructor() {
        this.renderer = new THREE.WebGLRenderer();
        this.renderer.setSize(300, 300);
        this.renderer.setClearColor(0x00ffff);
        this.renderer.clear();

        this.material = new THREE.RawShaderMaterial({
            uniforms: {
                time: { value: 1.0 }
            },
            vertexShader: `
            attribute vec3 position;

            void main() {
                gl_Position = vec4(position.xyz, 1.0);
            }
            `,
            fragmentShader: `
            void main() {
                gl_FragColor = vec4(0.5 + 0.5 * vec3(sin(gl_FragCoord.x), cos(gl_FragCoord.y), 0.0), 1.0);
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

const tracer = new WebGLTracer();
tracer.render();
document.body.append(tracer.renderer.domElement);