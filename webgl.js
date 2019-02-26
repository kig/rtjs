const mobile = /mobile/i.test(navigator.userAgent);
const dpr = (window.devicePixelRatio || 1);

class WebGLTracer {
    constructor(vgArray, traceGLSL) {
        var canvas = document.createElement( 'canvas' );
        var context = canvas.getContext( 'webgl2' );
        var texSize = Math.ceil(Math.sqrt(vgArray.length));
        console.log(texSize);
        if (texSize < 2048) {
            texSize = Math.pow(2, Math.ceil(Math.log2(texSize)));
        }
        console.log(texSize);
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
        camera.target = new THREE.Vector3(0, 1.01, 0);
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
                cameraInverseMatrix: { value: camera.inverseMatrix },
                deviceEpsilon: {value: mobile ? 0.01 : 0.0001},
                deviceEpsilonTrace: {value: mobile ? 0.05 : 0.01},
                costVis: {value: false},
                aaSize: {value: 1.0}
            },
            vertexShader: `
            #version 300 es

            precision highp float;
            precision highp int;

            void main() {
                gl_Position = vec4(position.xyz, 1.0);
            }
            `,
            fragmentShader: `
            #version 300 es

            precision highp float;
            precision highp int;
            
            ${traceGLSL}

            uniform float aaSize;

            out vec4 FragColor;    

            void main() {
                Array array = Array(arrayTexWidth);
                vec3 sum = vec3(0.0);
                for (float y = 0.0; y < aaSize; y++)
                for (float x = 0.0; x < aaSize; x++) {
                    vec3 c = trace(array, gl_FragCoord.xy + vec2(x,y) / aaSize);
                    c = -exp(-c * 0.5) + 1.0;
                    sum +=  c;
                }
                FragColor = vec4(sum / (aaSize*aaSize), 1.0);
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

        this.frame = 0;

        this.controls = new CameraControls(camera, canvas);
        camera.positionOffset.y = 0.1;
        this.startTime = Date.now();

        window.onresize = () => {
            this.controls.changed = true;
        }
    }

    render() {
        if (this.controls.changed) {
            this.controls.changed = false;

            const controlsActive = (this.controls.down || this.controls.pinching);

            this.material.uniforms.costVis.value = this.controls.debug;
            this.material.uniforms.aaSize.value = controlsActive ? 1 : 2;

            if (controlsActive) {
                this.renderer.setSize(window.innerWidth, window.innerHeight);
            } else {
                this.renderer.setSize(window.innerWidth*dpr, window.innerHeight*dpr);
            }
            
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
            if (this.frame === 0) {
                console.timeEnd('Load to first frame');
            }
            this.frame++;
        }
    }
}

(async function() {
    console.time('Load to first frame');

    const vgRes = await fetch('lib/voxelgrid.glsl');
    const traceRes = await fetch('lib/trace.glsl');
    const bunny = await ObjParse.load('bunny.obj');
    const vgText = await vgRes.text();
    const traceText = await traceRes.text();
        
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
    var grid = [4,8,4,4];
    if (bunnyTris.length < 10000) {
        // Use low-res grid
        // Fastest JS exec: [32]
        // Nice mix of VG steps + intersects: [4,4,4]
        // + Fast JS exec: [8, 8]
        grid = [32,4];
    }
    console.timeEnd('OBJ munging');

    console.time('Create VoxelGrid');
    const voxelGrid = new VoxelGrid3(bunnyTris.bbox.min, vec3(m), grid, 0);
    voxelGrid.addTriangles(bunnyTris);
    console.timeEnd('Create VoxelGrid');

    console.time('Serialize VoxelGrid');
    bunnyVG = voxelGrid.serialize();
    console.timeEnd('Serialize VoxelGrid');


    const tracer = new WebGLTracer(bunnyVG, vgText + '\n' + traceText);
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
