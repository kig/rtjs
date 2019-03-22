var CameraControls = function(camera, canvas) {
	var state = this;
	state.down = false;
	state.changed = false;
	state.wasDown = true;
	state.mouse = vec3(-1);
	state.alpha = Math.asin(0.5);
	state.theta = Math.PI*0.75;
	camera.positionOffset = new THREE.Vector3();
	this.camera = camera;

	var hammer = new Hammer(canvas);
	state.pinching = false;

	hammer.get('pinch').set({ enable: true });
	hammer.on('pinchstart', function(ev) {
		camera.startFov = camera.fov;
		state.pinching = true;
	});
	hammer.on('pinchend', function(ev) {
		state.pinching = false;
		state.changed = true;
		console.log('pinchend');
	});

	hammer.on('pinch', function(ev) {
		camera.fov = camera.startFov / ev.scale;
		camera.fov = max(min(camera.fov, 160), 1);
		state.pinching = true;
		state.changed = true;
		ev.preventDefault();
	});

	hammer.on('doubletap', function(ev) {
		state.debug = !state.debug;
		state.changed = true;
		ev.preventDefault();
	});

	canvas.onclick = function(ev) {
		if (!state.moved) {
			state.focusPoint = {x: ev.clientX, y: ev.clientY};
			state.changed = true;
		}
		ev.preventDefault();
	};

	canvas.onmousedown = function(ev) {
		if (ev.button === 0 || ev.button === 1) {
			state.moved = false;
			state.down = true;
			state.mouseStart = vec3(ev.clientX, ev.clientY, ev.button);
			state.mouse = vec3(ev.clientX, ev.clientY, ev.button);
			ev.preventDefault();
		}
	};
	canvas.onwheel = function(ev) {
		var wd = -ev.deltaY;
		wd = Math.min(Math.max(-15, wd), 15);
		var f = Math.pow(1.01, -wd);
		camera.fov *= f;
		camera.fov = max(min(camera.fov, 160), 1);
		state.changed = true;
		ev.preventDefault();
	};
	window.onmousemove = function(ev) {
		if (state.down) {
			var newMouse = vec3(ev.clientX, ev.clientY, state.mouse.button||0);
			var d = sub(newMouse, state.mouse);
			if (length(sub(state.mouseStart, newMouse)) > 5.0) {
				state.moved = true;
			}
			state.mouse = newMouse;
			if (ev.shiftKey || state.mouse.z > 0) {
				var mv = new THREE.Vector3(-d.x, 0, 0).multiplyScalar(0.005*camera.fov/55);
				camera.updateMatrixWorld();
				mv.applyMatrix3(camera.matrixWorld);
				mv.y += d.y * 0.005*camera.fov/55;
				camera.positionOffset.add(mv);
				camera.target.add(mv);
			} else {
				state.theta += d.x * 0.01;
				state.alpha += d.y * 0.01;
				state.alpha = max(0, min(state.alpha, Math.PI/2))
			}
			state.updateCameraPosition();
		}
	};
	window.onmouseup = function(ev) {
		if (state.down) {
			state.down = false;
			state.pinching = false;
			ev.preventDefault();
		}
	};
	canvas.addEventListener('touchstart', function(ev) {
		state.down = true;
		state.mouse = vec3(ev.touches[0].clientX, ev.touches[0].clientY, 0);
		ev.preventDefault();
	}, false);
	canvas.addEventListener('touchmove', function(ev) {
		if (state.down) {
			if (ev.touches.length > 1) {
				var newMouse = vec3(
					ev.touches.reduce((s,t) => s + t.clientX, 0) / ev.touches.length,
					ev.touches.reduce((s,t) => s + t.clientY, 0) / ev.touches.length,
					0
				);
				var d = sub(newMouse, state.mouse);
				state.mouse = newMouse;
				var mv = new THREE.Vector3(-d.x, 0, 0).multiplyScalar(0.01);
				camera.updateMatrixWorld();
				mv.applyMatrix3(camera.matrixWorld);
				mv.y += d.y * 0.01;
				camera.positionOffset.add(mv);
				camera.target.add(mv);
			} else {
				state.pinching = false;
				var newMouse = vec3(ev.touches[0].clientX, ev.touches[0].clientY, 0);
				var d = sub(newMouse, state.mouse);
				state.mouse = newMouse;
				state.theta += d.x * 0.01;
				state.alpha += d.y * 0.01;
				state.alpha = max(0, min(state.alpha, Math.PI/2));
			}
			state.updateCameraPosition();
			ev.preventDefault();
		}
	}, false);
	canvas.addEventListener('touchend', function(ev) {
		state.down = false;
		state.pinching = false;
		ev.preventDefault();
	}, false);
	canvas.addEventListener('touchcancel', function(ev) {
		state.down = false;
		state.pinching = false;
		ev.preventDefault();
	}, false);

	this.updateCameraPosition();
};

CameraControls.prototype.updateCameraPosition = function() {
	var state = this;
	state.camera.position.set(Math.cos(state.theta)*Math.cos(state.alpha), Math.sin(state.alpha), Math.cos(state.alpha)*Math.sin(state.theta)).normalize().multiplyScalar(4.5);
	state.camera.position.add(state.camera.positionOffset);
	state.changed = true;
};