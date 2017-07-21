var Ray = function(o, d, time) {
	this.o = o;
	this.d = d;
	this.transmit = vec3(1.0);
	this.light = vec3(0.0);
	this.time = time;
	this.bounce = 0;
};

var intersect = function(ray, scene) {
	var minT = 1/0;
	var obj = null;
	for (var i=0; i<scene.length; i++) {
		var t = scene[i].intersect(ray);
		if (t > 0 && t < minT) {
			minT = t;
			obj = scene[i];
		}
	}
	if (!obj) {
		return null;
	}
	return {distance: minT, obj: obj};
};

var setupRay = function(camera, x, y, w, h) {
	var uv = new THREE.Vector2(x/w*2 - 1, y/h*2 - 1);

	var origin = new THREE.Vector3();
	origin.setFromMatrixPosition( camera.matrixWorld );

	var direction = new THREE.Vector3();
	direction.set( uv.x, uv.y, 0.5 ).unproject( camera ).sub( origin ).normalize();

	var target = new THREE.Vector3();
	var focusDistance = origin.distanceTo(camera.focusPoint);
	target.addVectors( origin, direction.multiplyScalar(focusDistance));

	origin.add( new THREE.Vector3().copy(diskPoint()).multiplyScalar(camera.apertureSize).applyMatrix4(camera.matrixWorld) );
	direction.subVectors( target, origin ).normalize();

	var ray = new Ray(vec3(origin.x, origin.y, origin.z), vec3(direction.x, direction.y, direction.z), 0);
	return ray;
};


var loader = new THREE.OBJLoader;
loader.load('bunny.obj', function(bunny) {
	var camera = new THREE.PerspectiveCamera(55, 1, 0.1, 100);
	camera.position.set(3.5, 1.5, 3.5).normalize().multiplyScalar(4.0);
	camera.target = new THREE.Vector3(0, 0.5, 0);
	camera.lookAt(camera.target);
	camera.updateProjectionMatrix();
	camera.updateMatrixWorld();

	var canvas = document.createElement('canvas');
	document.body.appendChild(canvas);
	var ctx = canvas.getContext('2d');

	var scene = [
		new Sphere(vec3(2.5,1,1), 1, vec3(1.0, 0.7, 0.3))
	];

	bunny.children[0].geometry.computeBoundingBox();
	var bbox = bunny.children[0].geometry.boundingBox;
	var verts = bunny.children[0].geometry.attributes.position.array;
	var normals = bunny.children[0].geometry.attributes.normal.array;

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

	var bunnyTris = [];

	for (var i = 0; i < verts.length; i += 3*3) {
		var u = vec3(verts[i], verts[i+1], verts[i+2]);
		var v = vec3(verts[i+3], verts[i+4], verts[i+5]);
		var w = vec3(verts[i+6], verts[i+7], verts[i+8]);
		var x = vec3(normals[i], normals[i+1], normals[i+2]);
		var y = vec3(normals[i+3], normals[i+4], normals[i+5]);
		var z = vec3(normals[i+6], normals[i+7], normals[i+8]);
		bunnyTris.push(new Triangle([u,v,w], [x,y,z], add(vec3(0.8), mulS(absV(u.xxy),0.2))));
	}

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

	var acceleration = "BVH";
	var paused = false;

	window.useBVH.onclick = function() {
		acceleration = "BVH";
		wasDown = true;
	};

	window.useVoxelGrid.onclick = function() {
		acceleration = "VoxelGrid";
		wasDown = true;
	};

	window.pauseButton.onclick = function() {
		paused = !paused;
	};

	window.canvasSize.onchange = window.aaSize.onchange = window.apertureSize.onchange = function() {
		wasDown = true;
	};

	var down = false;
	var wasDown = true;
	var mouse = vec3(-1);
	var alpha = Math.asin(0.5);
	var theta = 0;
	canvas.onmousedown = function(ev) {
		if (ev.button === 0) {
			down = true;
			mouse = vec3(ev.clientX, ev.clientY, ev.button);
			ev.preventDefault();
		}
	};
	window.onmousemove = function(ev) {
		if (down) {
			console.log(ev);
			var newMouse = vec3(ev.clientX, ev.clientY, ev.button||0);
			var d = sub(newMouse, mouse);
			mouse = newMouse;
			theta += d.x * 0.01;
			alpha += d.y * 0.01;
			alpha = max(0, min(alpha, Math.PI/2))
		}
	};
	window.onmouseup = function(ev) {
		if (down) {
			down = false;
			ev.preventDefault();
		}
	};
	canvas.addEventListener('touchstart', function(ev) {
		down = true;
		mouse = vec3(ev.touches[0].clientX, ev.touches[0].clientY, 0);
		ev.preventDefault();
	}, false);
	canvas.addEventListener('touchmove', function(ev) {
		if (down) {
			var newMouse = vec3(ev.touches[0].clientX, ev.touches[0].clientY, 0);
			var d = sub(newMouse, mouse);
			mouse = newMouse;
			theta += d.x * 0.01;
			alpha += d.y * 0.01;
			alpha = max(0, min(alpha, Math.PI/2))
			ev.preventDefault();
		}
	}, false);
	canvas.addEventListener('touchend', function(ev) {
		down = false;
		ev.preventDefault();
	}, false);
	canvas.addEventListener('touchcancel', function(ev) {
		down = false;
		ev.preventDefault();
	}, false);

	var t = 0;

	function render() {
		if (!down && !wasDown && paused) return;
		if (!paused) t += 16;
		camera.position.set(Math.cos(theta)*Math.cos(alpha), Math.sin(alpha), Math.cos(alpha)*Math.sin(theta)).normalize().multiplyScalar(4.5);
		camera.target = new THREE.Vector3(0, 0.75, 0);
		camera.focusPoint = vec3(-1, 0.7, 0.2);
		camera.lookAt(camera.target);
		camera.updateProjectionMatrix();
		camera.updateMatrixWorld();

		var canvasSize = parseInt(window.canvasSize.value);
		var AA_SIZE = parseInt(window.aaSize.value);
		if (down) {
			canvasSize = 50;
			AA_SIZE = 1;
			wasDown = true;
		} else if (wasDown) {
			wasDown = false;
		}
		var apertureSize = parseInt(window.apertureSize.value);
		camera.apertureSize = Math.pow(1.33, -apertureSize);

		canvas.width = canvas.height = canvasSize;

		var id = ctx.getImageData(0, 0, canvasSize, canvasSize);

		window.debug.innerHTML = "";

		scene.forEach(function(o,i) {
		 	o.center = vec3(o.center.x, abs(sin(i + t/300))*o.radius*2 + o.radius, o.center.z);
		 });

		if (acceleration === 'VoxelGrid') {

			(console || window.console).time("voxelGrid build");

			var voxelGrid = new VoxelGrid2(8, vec3(-1.1,-0.1,-1.1), vec3(2.2), 1, 8);
			bunnyTris.forEach(function(o) {
				voxelGrid.add(o);
			});

			(console || window.console).timeEnd("voxelGrid build");

			(console || window.console).time("voxelGrid2 build");

			var voxelGrid2 = new VoxelGrid2(16, vec3(-8.1,-0.1,-8.1), vec3(16.2), 0, 1);
			scene.forEach(function(o,i) {
				voxelGrid2.add(o);
			});

			(console || window.console).timeEnd("voxelGrid2 build");

		} else if (acceleration === 'BVH') {

			(console || window.console).time("BVH build");

			var bvh = new BVHNode(bunnyTris.concat(scene));
			console.log("Built BVH, size", bvh.subTreeSize);

			(console || window.console).timeEnd("BVH build");

		}

		console.time("trace");

		var rays = [];
		var epsilon = 0.0001;

		for (var y=0; y<canvasSize; y++) {
			for (var x=0; x<canvasSize; x++) {
				for (var dy=0; dy<AA_SIZE; dy++) {
					for (var dx=0; dx<AA_SIZE; dx++) {
						var r = setupRay(camera, x+dx/AA_SIZE, (canvasSize-y-1)+dy/AA_SIZE, canvasSize, canvasSize);
						rays.push(r);
					}
				}
			}
		}

		var plane = new Plane(vec3(0,0,0), vec3(0,1,0), vec3(0.5));
		var rayCount = 0;
		var lastRayCount = rayCount;
		VoxelGrid.stepCount = 0;
		VoxelGrid.cmpCount = 0;
		BVHNode.visitedCount = 0;
		BVHNode.primitiveTests = 0;
		console.log("Tracing " + rays.length + " primary rays");
		for (var j=0; j<6; j++) {
			for (var i=0; i<rays.length; i++) {
				var r = rays[i];
				if (r.finished) continue;
				rayCount++;
				if (acceleration === 'BVH') {
					var hit = bvh.intersect(r);	
				} else {
					var hit = voxelGrid.intersect(r);	
					var hit1 = voxelGrid2.intersect(r);
					if (!hit || (hit1 && hit1.distance < hit.distance)) {
						hit = hit1;
					}
				}
				var hit2 = plane.intersect(r);
				if (hit2 > 0 && (!hit || hit2 < hit.distance)) {
					hit = {distance: hit2, obj: plane};
				}
				if (hit) {
					if (r.bounce < 5) {
						r.o = add(r.o, mulS(r.d, hit.distance));
						var c = hit.obj.color(r);
						r.transmit = mul(r.transmit, c);
						var nml = hit.obj.normal(r.o);
						r.d = normalize(reflect(r.d, nml));
						r.o = add(r.o, mulS(r.d, epsilon));
						r.bounce++;
					} else {
						r.finished = true;
					}
				} else {
					var bg = mulS(vec3(0.4+sat(-r.d.y), 0.6, 0.8+sat(r.d.y)*abs(r.d.z)), 2);
					bg = add(bg, mulS(vec3(10.0, 6.0, 4.0), pow(sat(dot(r.d, normalize(vec3(5.0, 15.0, 10.0)))), 64.0) ));
					r.light = add(r.light, mul(r.transmit, bg));
					r.finished = true;
				}
			}
			console.log("Bounce " + j + ", traced " + (rayCount-lastRayCount));
			lastRayCount = rayCount;
		}

		console.timeEnd("trace");

		console.log("Traced " + rayCount + " rays");

		
		if (acceleration === 'VoxelGrid') {
			console.log(VoxelGrid.stepCount, "VoxelGrid steps");
			console.log(VoxelGrid.cmpCount, "VoxelGrid primitive intersection tests");
			console.log(VoxelGrid.stepCount / rayCount, "VG steps per ray")
			console.log(VoxelGrid.cmpCount / rayCount, "VG primitive intersection tests per ray");
		}

		if (acceleration === 'BVH') {
			console.log(BVHNode.visitedCount, "BVH visited nodes");
			console.log(BVHNode.primitiveTests, "BVH primitive intersection tests");
			console.log(BVHNode.visitedCount / rayCount, "BVH nodes per ray")
			console.log(BVHNode.primitiveTests / rayCount, "BVH primitive intersection tests per ray");
		}

		for (var i=0; i<canvasSize*canvasSize; i++) {
			var c = vec3();
			for (var dy=0; dy<AA_SIZE; dy++) {
				for (var dx=0; dx<AA_SIZE; dx++) {
					var r = rays[i*AA_SIZE*AA_SIZE+dy*AA_SIZE+dx];
					c = add(c, addS(neg(expV(neg(r.light))), 1));
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

