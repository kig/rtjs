
// Two-level binary voxel grid to accelerate triangle search

uniform mediump isampler2D voxelIndex;
uniform mediump usampler2D triIndices;
uniform sampler2D triangles;
uniform sampler2D normals;

uniform highp int voxelIndexWidth;
uniform highp int triIndicesWidth;
uniform highp int trianglesWidth;
uniform highp int normalsWidth;

uniform vec3 vgOrigin;
uniform float vgScale;
uniform highp int vgSize;

void intersectTris(inout Ray ray, in int coff, inout Hit closestHit) {
    for (int j = 0; j < 16; j++) {
        uint triIndex = readUint(triIndices, triIndicesWidth, coff + j);
        if (triIndex == 0u) {
            break;
        }
        if (ray.lastTested != int(triIndex)) {
            if (costVis) {
                ray.light.r += 0.1;
            }
            intersectTri(triangles, trianglesWidth, ray, int(triIndex - 1u), closestHit);
        }
    }
}

void getSubNode(
    in Ray ray, in ivec3 cstep, in vec3 deltaDist, in int size, 
    inout ivec3 c, inout vec3 next, inout int ci, 
    inout int voxelIndexOffset, inout int voxelsOffset, 
    inout vec3 origin, inout float scale
) {
    voxelIndexOffset = voxelIndexOffset + (((1+ci) * size*size*size) >> 4);
    voxelsOffset = voxelsOffset + (ci * size*size*size) * 16;

    origin = origin + (vec3(c) * scale);

    float t = forceIntersectBox(ray, origin, vec3(scale));

    scale = scale / float(size);

	vec3 cf = toGrid(ray.o + ray.d * (t + scale * 0.001), scale, origin);
	c = ivec3(cf);

	ci = c.z * size * size + c.y * size + c.x;

    next = (max(vec3(0.0), vec3(cstep)) + vec3(cstep) * (vec3(c) - cf)) * deltaDist;
}

void intersectGrid(inout Ray ray, inout Hit closestHit) {
    vec3 origin = vgOrigin;
    int size = vgSize;
    float scale = vgScale;

    int voxelIndexOffset = 0;
    int voxelsOffset = 0;

    float t = intersectBox(ray, origin, vec3(float(size) * scale));
    if (t < 0.0) {
        return;
    }
    
    vec3 cp = ray.o + ray.d * (t + 0.001 * scale);
    vec3 cf = toGrid(cp, scale, origin);
    ivec3 c = ivec3(cf);

    vec3 deltaDist = vec3(
        length(ray.d * ray.invD.x),
        length(ray.d * ray.invD.y),
        length(ray.d * ray.invD.z)
    );
    ivec3 cstep = ivec3(sign(ray.d));
    vec3 next = (max(vec3(0.0), vec3(cstep)) + vec3(cstep) * (vec3(c) - cf)) * deltaDist;

    int ci = c.z * size * size + c.y * size + c.x;

    ivec3 topC = c;
    vec3 topNext = next;
    int topCi = ci;
    vec3 topOrigin = origin;
    float topScale = scale;

    bool onTopLevel = true;

    // Step through the grid while we're inside it
    for (int i = 0; i < 3*size*size; i++) {
        if (c.x < 0 || c.y < 0 || c.z < 0 || c.x >= size || c.y >= size || c.z >= size) {
            if (onTopLevel) {
                // escaped top level
                break;
            } else {
                // pop stack
                onTopLevel = true;
                c = topC;
                next = topNext;
                ci = topCi;
                origin = topOrigin;
                scale = topScale;
                voxelIndexOffset = 0;
                voxelsOffset = 0;
            }
        } else {
            if (costVis) {
                ray.light.b += 0.01;
            }
            int vi = (readInt(voxelIndex, voxelIndexWidth, voxelIndexOffset + (ci >> 4)) >> (ci & 15)) & 1;
            if (vi != 0) {
                // When at top-level node, jump to second-level node and continue intersect
                // (push current node on stack, replace with sub-node stuff)
                if (onTopLevel) {
                    topC = c;
                    topNext = next;
                    topCi = ci;
                    getSubNode(ray, cstep, deltaDist, size, c, next, ci, voxelIndexOffset, voxelsOffset, origin, scale);
                    onTopLevel = false;
                    continue;
                } else {
                    intersectTris(ray, voxelsOffset + ci * 16, closestHit);
                    if (closestHit.index >= 0) {
                        ivec3 p = ivec3(toGrid(ray.o + ray.d * closestHit.distance, scale, origin));
                        if (all(equal(p, c))) {
                            return;
                        }
                    }
                }
            }
        }
        if (next.x < next.y) {
            if (next.x < next.z) {
                next.x += deltaDist.x;
                c.x += cstep.x;
                ci += cstep.x;
            } else {
                next.z += deltaDist.z;
                c.z += cstep.z;
                ci += cstep.z * size * size;
            }
        } else if (next.y < next.z) {
            next.y += deltaDist.y;
            c.y += cstep.y;
            ci += cstep.y * size;
        } else {
            next.z += deltaDist.z;
            c.z += cstep.z;
            ci += cstep.z * size * size;
        }
    }
}
