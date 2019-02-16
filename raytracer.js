var Ray = function(o, d, time) {
	this.o = o;
	this.d = d;
	this.transmit = vec3(1.0);
	this.light = vec3(0.0);
	this.time = time;
	this.pathLength = 0;
	this.bounce = 0;
	this.finished = false;
};

var intersect = function(ray, scene) {
	var minT = 1/0;
	var hit = null;
	for (var i=0; i<scene.length; i++) {
		var t = scene[i].intersect(ray);
		if (t && t.distance < minT) {
			minT = t.distance;
			hit = t;
		}
	}
	return hit;
};

const origin = new THREE.Vector3();
const direction = new THREE.Vector3();
const target = new THREE.Vector3();
const diskPointV = new THREE.Vector3();
const uv = new THREE.Vector2();
var setupRay = function(ray, camera, x, y, w, h) {
	uv.set(x/w*2 - 1, y/h*2 - 1);

	origin.setFromMatrixPosition( camera.matrixWorld );

	direction.set( uv.x, uv.y, 0.5 ).unproject( camera ).sub( origin ).normalize();

	const focusDistance = origin.distanceTo(camera.focusPoint);
	target.addVectors( origin, direction.multiplyScalar(focusDistance));

	origin.add( diskPointInPlace(diskPointV).multiplyScalar(camera.apertureSize).applyMatrix4(camera.matrixWorld) );
	direction.subVectors( target, origin ).normalize();

	ray.o.x = origin.x;
	ray.o.y = origin.y;
	ray.o.z = origin.z;
	ray.d.x = direction.x;
	ray.d.y = direction.y;
	ray.d.z = direction.z;
	ray.transmit.x = ray.transmit.y = ray.transmit.z = 1;
	ray.light.x = ray.light.y = ray.light.z = 0;
	ray.pathLength = 0;
	ray.bounce = 0;
	ray.finished = false;
};

var trace = function(rays, raysLength, scene, console) {
	const epsilon = 0.0001;
	console.time("trace");

	console.log("Tracing " + raysLength + " primary rays");
	const plane = new Plane(vec3(0,0,0), vec3(0,1,0), vec3(0.5));
	let rayCount = 0;
	let lastRayCount = rayCount;
	for (let j=0; j<6; j++) {
		for (let i=0; i<raysLength; i++) {
			const r = rays[i];
			if (r.finished) continue;
			rayCount++;
			let hit = scene.intersect(r);
			const hit2 = plane.intersect(r);
			if (hit2 && (!hit || hit2.distance < hit.distance)) {
				hit = hit2;
			}
			if (hit) {
				r.o = add(r.o, mulS(r.d, hit.distance));
				r.pathLength += hit.distance;
				const c = hit.obj.color(r);
				r.transmit = mul(r.transmit, c);
				const nml = hit.obj.normal(r.o);
				r.d = normalize(add(reflect(r.d, nml), mulS(vec3(Math.random()-.5, Math.random()-.5, 2*(Math.random()-.5)), 0.1)));
				r.o = add(r.o, mulS(nml, epsilon));
				r.bounce++;
			} else {
				let bg = mulS(vec3(0.6+sat(-r.d.y), 0.7, 0.8+(0.4*r.d.x)*abs(r.d.z)), 1);
				bg = add(bg, mulS(vec3(10.0, 6.0, 4.0), 4*pow(sat(dot(r.d, normalize(vec3(6.0, 10.0, 8.0)))), 64.0) ));
				bg = add(bg, mulS(vec3(3, 5, 7), abs(1-r.d.z)));
				r.light = mix(add(r.light, mul(r.transmit, bg)), bg, 1-Math.exp(-r.pathLength/40));
				r.finished = true;
			}
		}
		console.log("Bounce " + j + ", traced " + (rayCount-lastRayCount));
		lastRayCount = rayCount;
	}
	for (let i=0; i<raysLength; i++) {
		const r = rays[i];
		if (r.finished) continue;
		var bg = mulS(vec3(0.6+sat(-r.d.y), 0.7, 0.8+(0.4*r.d.x)*abs(r.d.z)), 1);
		bg = add(bg, mulS(vec3(10.0, 6.0, 4.0), 4*pow(sat(dot(r.d, normalize(vec3(6.0, 10.0, 8.0)))), 64.0) ));
		bg = add(bg, mulS(vec3(3, 5, 7), abs(1-r.d.z)));
		r.light = mulS(bg, 1-Math.exp(-r.pathLength/40));
		r.finished = true;
	}

	console.timeEnd("trace");

	console.log("Traced " + rayCount + " rays");
	return rayCount;
};

var cachedAcceleration = {};
var getAcceleration = function(bunnyTris, scene, bvhWidth, acceleration, rays, cache, console) {
	if (cache && cachedAcceleration[acceleration]) {
		return cachedAcceleration[acceleration];
	}

	if (acceleration === 'PathCache') {
		console.log("PathCache is a ray bounce processing baseline. It traces the scene using a BVH and stores the results. On the PathCache trace, it replays the results. Think of it as the optimal acceleration structure, getting the nearest intersecting object for a ray in O(1).");
		console.time("PathCache build");
		var bvh = new BVHNode(bunnyTris.concat(scene), bvhWidth);

		trace(rays, bvh, {time:function(){}, timeEnd:function() {}, log:function(){}}, function(r, hit) {
			if (!r.hits) {
				r.hitIndex = 0;
				r.hits = [];
				r.oo = addS(r.o, 0);
				r.od = addS(r.d, 0);
			}
			r.hits.push(hit);
		});
		rays.forEach(function(r) {
			r.o = r.oo;
			r.d = r.od;
			r.light = vec3(0);
			r.transmit = vec3(1);
			r.finished = false;
			r.bounce = 0;
		});
		var accel = {
			intersect: function(r) {
				return r.hits[r.hitIndex++];
			}
		};
		console.timeEnd("PathCache build");
	}

	if (acceleration === 'VoxelGrid') {

		console.time("voxelGrid build");

		var size = sub(bunnyTris.bbox.max, bunnyTris.bbox.min);
		var m = Math.max(size.x, size.y, size.z);
		var grid = [8,4,4,4];
		if (bunnyTris.length < 10000) {
			// Use low-res grid
			// Fastest JS exec: [64]
			// Nice mix of VG steps + intersects: [4,4,4]
			// + Fast JS exec: [8, 8]
			grid = [8, 8];
		}
		const voxelGrid = new VoxelGrid3(bunnyTris.bbox.min, vec3(m), grid, 0);
		voxelGrid.addTriangles(bunnyTris);
		// for (var i = 0; i < bunnyTris.length; i++) {
		// 	voxelGrid.add(bunnyTris[i]);
		// }

		const blob = voxelGrid.serialize();
		window.console.log(blob);
		
		// var accel = voxelGrid;
		var accel = new SerializedVG(blob, bunnyTris[0]._color);

		console.timeEnd("voxelGrid build");

		// console.time("voxelGrid2 build");

		// var voxelGrid2 = new VoxelGrid2(16, vec3(-8.1,-0.1,-8.1), vec3(16.2), 0, [1], 0);
		// scene.forEach(function(o,i) {
		// 	voxelGrid2.add(o);
		// });

		// var accel = voxelGrid;
		// var x = {
		// 	intersect: function(r) {
		// 		var hit = voxelGrid.intersect(r);	
		// 		var hit1 = voxelGrid2.intersect(r);
		// 		if (!hit || (hit1 && hit1.distance < hit.distance)) {
		// 			hit = hit1;
		// 		}
		// 		return hit;
		// 	}
		// };

		// console.timeEnd("voxelGrid2 build");

	} else if (acceleration === 'BeamSphere') {

		window.console.time("BeamSphere build");

		window.console.time("BeamSphere init");
		var size = sub(bunnyTris.bbox.max, bunnyTris.bbox.min);
		var mid = add(bunnyTris.bbox.min, mulS(size, 0.5));
		var radius = length(size) / 2;
		var accel = new BeamSphere(mid, radius, 5);
		window.console.timeEnd("BeamSphere init");

		for (var i = 0; i < bunnyTris.length; i++) {
			accel.add(bunnyTris[i]);
		}
		accel.sort();

		window.console.timeEnd("BeamSphere build");
		window.console.log("Added triangles", BeamSphere.addedTriangles);

	} else if (acceleration === 'BVH') {

		console.time("BVH build");

		var accel = new BVHNode(bunnyTris, bvhWidth);
		// var bvh = new BVHNode(bunnyTris, bvhWidth);
		// var objs = [];
		// for (var i=0; i<1; i++) {
		// 	var tn = new TransformNode(bvh);
		// 	tn.transform.position.copy(randomVec3Unit());
		// 	tn.transform.rotation.y = random() * 2 * Math.PI;
		// 	tn.transform.updateMatrix();
		// 	tn.transform.inverseMatrix = new THREE.Matrix4();
		// 	tn.transform.inverseMatrix.getInverse(tn.transform.matrix);
		// 	objs.push(tn);
		// }
		// var accel = new BVHNode(objs, 2);
		// accel = bvh;

		console.log("Built BVH, width", bvhWidth, "size", accel.subTreeSize);

		console.timeEnd("BVH build");

	}
	if (cache) {
		cachedAcceleration[acceleration] = accel;
	}
	return accel;
};

ObjParse.load('bunny.obj').then(function(bunny) {
	var camera = new THREE.PerspectiveCamera(55, 1, 0.1, 100);
	camera.target = new THREE.Vector3(0, 0.75, 0);
	camera.focusPoint = vec3(-1, 0.7, 0.2);
	camera.positionOffset = new THREE.Vector3();
	camera.lookAt(camera.target);
	camera.updateProjectionMatrix();
	camera.updateMatrixWorld();

	var canvas = document.createElement('canvas');
	document.body.appendChild(canvas);
	var ctx = canvas.getContext('2d');

	var scene = [
		new Sphere(vec3(2.5,1,1), 1, vec3(1.0, 0.7, 0.3))
	];

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

	for (var i=0; i<100; i++) {
		var c = mulS(randomVec3Unit(), 4 + random() * 3);
		var r = 0.05 + 0.2 * random();
		c.y = r;
		var color = randomVec3Positive();
		scene.push(new Sphere(c, r, color));
	}

	var console = {
		timers: {},
		time: function(n) {
			this.timers[n] = performance.now();
		},
		timeEnd: function(n) {
			var ms = performance.now() - this.timers[n];
			this.log(n + " " + ms + " ms");
		},
		log: function() {
			var d = document.createElement('div');
			d.textContent = [].join.call(arguments, " ");
			window.debug.appendChild(d);
		}
	};

	var acceleration = "VoxelGrid";
	var paused = true;

	window.useBVH.onclick = function() {
		acceleration = "BVH";
		controls.wasDown = true;
	};

	window.useVoxelGrid.onclick = function() {
		acceleration = "VoxelGrid";
		controls.wasDown = true;
	};

	window.useBeamSphere.onclick = function() {
		acceleration = "BeamSphere";
		controls.wasDown = true;
	};
	
	window.usePathCache.onclick = function() {
		acceleration = "PathCache";
		controls.wasDown = true;
	};
	
	window.pauseButton.onclick = function() {
		paused = !paused;
	};

	window.canvasSize.onchange = 
	window.aaSize.onchange = 
	window.apertureSize.onchange =
	window.bvhWidth.onchange = 
	function() {
		controls.wasDown = true;
	};

	var controls = new CameraControls(camera, canvas);


	var t = 0;

	var rays = [];

	function render() {
		if (!controls.changed && !controls.down && !controls.wasDown && paused) return;
		if (!paused) t += 16;
		camera.lookAt(camera.target);
		camera.updateProjectionMatrix();
		camera.updateMatrixWorld();

		var canvasSize = parseInt(window.canvasSize.value);
		var AA_SIZE = parseInt(window.aaSize.value);
		if (controls.down || controls.changed) {
			canvasSize = 50;
			AA_SIZE = 1;
			controls.wasDown = true;
			controls.changed = false;
		} else if (controls.wasDown) {
			controls.wasDown = false;
		}
		var apertureSize = parseInt(window.apertureSize.value);
		camera.apertureSize = Math.pow(1.33, -apertureSize);
		if (apertureSize === 10) {
			camera.apertureSize = 0;
		}

		var bvhWidth = parseInt(window.bvhWidth.value);

		canvas.width = canvas.height = canvasSize;

		var id = ctx.getImageData(0, 0, canvasSize, canvasSize);

		window.debug.innerHTML = "";


		const totalRayCount = canvasSize*canvasSize*AA_SIZE*AA_SIZE;
		if (rays.length < totalRayCount) {
			rays = [];
			for (var i=0; i<totalRayCount; i++) {
				rays[i] = new Ray(vec3(0), vec3(0), 0);
			}
		}

		(console || window.console).time("Create rays");

		for (var y=0, i=0; y<canvasSize; y++) {
			for (var x=0; x<canvasSize; x++) {
				for (var dy=0; dy<AA_SIZE; dy++) {
					for (var dx=0; dx<AA_SIZE; dx++, i++) {
						setupRay(rays[i], camera, x+dx/AA_SIZE, (canvasSize-y-1)+dy/AA_SIZE, canvasSize, canvasSize);
					}
				}
			}
		}

		(console || window.console).timeEnd("Create rays");


		scene.forEach(function(o,i) {
		 	o.center = vec3(o.center.x, abs(sin(i + t/300))*o.radius*2 + o.radius, o.center.z);
		});

		var accel = getAcceleration(bunnyTris, scene, bvhWidth, acceleration, rays, true, console);


		VoxelGrid.stepCount = 0;
		VoxelGrid.cmpCount = 0;
		BeamSphere.sphereTests = 0;
		BeamSphere.beamTests = 0;
		BeamSphere.primitiveTests = 0;
		BVHNode.visitedCount = 0;
		BVHNode.primitiveTests = 0;

		var rayCount = trace(rays, totalRayCount, accel, console);
		
		if (acceleration === 'VoxelGrid') {
			console.log(VoxelGrid.stepCount, "VoxelGrid steps");
			console.log(VoxelGrid.cmpCount, "VoxelGrid primitive intersection tests");
			console.log(VoxelGrid.stepCount / rayCount, "VG steps per ray");
			console.log(VoxelGrid.cmpCount / rayCount, "VG primitive intersection tests per ray");
		}

		if (acceleration === 'BeamSphere') {
			console.log(BeamSphere.sphereTests, "BeamSphere sphere tests");
			console.log(BeamSphere.beamTests, "BeamSphere beam tests");
			console.log(BeamSphere.primitiveTests, "BeamSphere primitive intersection tests");
			console.log(BeamSphere.sphereTests / rayCount, "BS sphere tests per ray");
			console.log(BeamSphere.beamTests / rayCount, "BS beam tests per ray");
			console.log(BeamSphere.primitiveTests / rayCount, "BS primitive intersection tests per ray");
		}

		if (acceleration === 'BVH') {
			console.log(BVHNode.visitedCount, "BVH visited nodes");
			console.log(BVHNode.primitiveTests, "BVH primitive intersection tests");
			console.log(BVHNode.visitedCount / rayCount, "BVH nodes per ray");
			console.log(BVHNode.primitiveTests / rayCount, "BVH primitive intersection tests per ray");
		}

		for (var i=0; i<canvasSize*canvasSize; i++) {
			var c = vec3();
			for (var dy=0; dy<AA_SIZE; dy++) {
				for (var dx=0; dx<AA_SIZE; dx++) {
					var r = rays[i*AA_SIZE*AA_SIZE+dy*AA_SIZE+dx];
					c = add(c, addS(neg(expV(neg(mulS(r.light, 0.5)))), 1));
				}
			}
			c = mulS(c, 1/(AA_SIZE*AA_SIZE));
			id.data[i*4 + 0] = c.x * 255;
			id.data[i*4 + 1] = c.y * 255;
			id.data[i*4 + 2] = c.z * 255;
			id.data[i*4 + 3] = 255;
		}

		ctx.putImageData(id, 0, 0);
	};

	var tick = function() {
		render();
		requestAnimationFrame(tick);
	};

	requestAnimationFrame(tick);
});
