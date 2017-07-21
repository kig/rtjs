var CameraControls = function(camera, canvas) {
	var state = this;
	state.down = false;
	state.wasDown = true;
	state.mouse = vec3(-1);
	state.alpha = Math.asin(0.5);
	state.theta = 0;
	camera.positionOffset = new THREE.Vector3();
	this.camera = camera;

	canvas.onmousedown = function(ev) {
		if (ev.button === 0 || ev.button === 1) {
			state.down = true;
			state.mouse = vec3(ev.clientX, ev.clientY, ev.button);
			ev.preventDefault();
		}
	};
	canvas.onmousewheel = function(ev) {
		var wd = ev.wheelDeltaY;
		var f = Math.pow(1.01, -wd);
		camera.fov *= f;
		camera.fov = max(min(camera.fov, 160), 5);
		state.wasDown = true;
	};
	window.onmousemove = function(ev) {
		if (state.down) {
			var newMouse = vec3(ev.clientX, ev.clientY, state.mouse.button||0);
			var d = sub(newMouse, state.mouse);
			state.mouse = newMouse;
			if (ev.shiftKey || state.mouse.z > 0) {
				var mv = new THREE.Vector3(-d.x, 0, 0).multiplyScalar(0.01);
				camera.updateMatrixWorld();
				mv.applyMatrix3(camera.matrixWorld);
				mv.y += d.y * 0.01;
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
			var newMouse = vec3(ev.touches[0].clientX, ev.touches[0].clientY, 0);
			var d = sub(newMouse, state.mouse);
			state.mouse = newMouse;
			state.theta += d.x * 0.01;
			state.alpha += d.y * 0.01;
			state.alpha = max(0, min(state.alpha, Math.PI/2));
			state.updateCameraPosition();
			ev.preventDefault();
		}
	}, false);
	canvas.addEventListener('touchend', function(ev) {
		state.down = false;
		ev.preventDefault();
	}, false);
	canvas.addEventListener('touchcancel', function(ev) {
		state.down = false;
		ev.preventDefault();
	}, false);

	this.updateCameraPosition();
};

CameraControls.prototype.updateCameraPosition = function() {
	var state = this;
	state.camera.position.set(Math.cos(state.theta)*Math.cos(state.alpha), Math.sin(state.alpha), Math.cos(state.alpha)*Math.sin(state.theta)).normalize().multiplyScalar(4.5);
	state.camera.position.add(state.camera.positionOffset);
};