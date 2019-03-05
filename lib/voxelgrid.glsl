
// Two-level voxel grid to accelerate triangle search


vec3 fromGrid(vec3 coord, float scale, vec3 origin) {
    return origin + (coord * scale);
}

vec3 toGrid(vec3 point, float scale, vec3 origin) {
    return (point - origin) / scale;
}

void intersectTris(in Array array, inout Ray ray, in int coff, in int childSize, inout Hit closestHit) {
    for (int j = 0; j < childSize; j++) {
        int triIndex = readInt(array, coff + j);
        if (triIndex < 0) {
            break;
        }
        if (ray.lastTested != triIndex) {
            if (costVis) {
                ray.light.r += 0.1;
            }
            intersectTri(array, ray, triIndex, closestHit);
        }
    }
}

#ifndef BREADTH_FIRST

void intersectGridLeaf(in Array array, inout Ray ray, in int headOff, inout Hit closestHit, 
                        ivec3 cstep, vec3 deltaDist, vec3 origin, int size, float scale, vec3 cp) {
    // vec3 origin = readVec3(array, headOff);
    // int size = readInt(array, headOff + 3);
    // vec3 dims = readVec3(array, headOff + 4);
    int childSize = readInt(array, headOff + 7);

    int voxelsOff = headOff + 8;
    int childIndexOff = voxelsOff + size * size * size;

    vec3 cf = toGrid(cp, scale, origin);
    ivec3 c = ivec3(cf);

    if (c.x < 0 || c.y < 0 || c.z < 0 || c.x >= size || c.y >= size || c.z >= size) {
        return;
    }

    vec3 next = (max(vec3(0.0), vec3(cstep)) + vec3(cstep) * (vec3(c) - cf)) * deltaDist;

    int ci = c.z * size * size + c.y * size + c.x;

    // Step through the grid while we're inside it
    for (int i = 0; i < 3*size; i++) {
        if (costVis) {
            ray.light.b += 0.01;
        }
        int vi = readInt(array, voxelsOff + ci);
        if (vi > 0) {
            if (costVis) {
                ray.light.g += 0.1;
            }
            int coff = childIndexOff + (vi - 1) * childSize;
            intersectTris(array, ray, coff, childSize, closestHit);
            if (closestHit.index >= 0) {
                ivec3 p = ivec3(toGrid(ray.o + ray.d * closestHit.distance, scale, origin));
                if (all(equal(p, c))) {
                    return;
                }
            }
        }
        if (next.x < next.y) {
            if (next.x < next.z) {
                next.x += deltaDist.x;
                c.x += cstep.x;
                ci += cstep.x;
                if (c.x < 0 || c.x >= size) {
                    return;
                }
            } else {
                next.z += deltaDist.z;
                c.z += cstep.z;
                ci += cstep.z * size * size;
                if (c.z < 0 || c.z >= size) {
                    return;
                }
            }
        } else if (next.y < next.z) {
            next.y += deltaDist.y;
            c.y += cstep.y;
            ci += cstep.y * size;
            if (c.y < 0 || c.y >= size) {
                return;
            }
        } else {
            next.z += deltaDist.z;
            c.z += cstep.z;
            ci += cstep.z * size * size;
            if (c.z < 0 || c.z >= size) {
                return;
            }
        }
    }
}

void intersectGridNode(in Array array, inout Ray ray, in int headOff, inout Hit closestHit) {
    vec3 origin = readVec3(array, headOff);
    int osize = readInt(array, headOff + 3);
    int size = abs(osize);
    bool isLeaf = osize > 0;
    vec3 dims = readVec3(array, headOff + 4);
    int childSize = readInt(array, headOff + 7);
    float scale = dims.x / float(size);
    int voxelsOff = headOff + 8;
    int childIndexOff = voxelsOff + size * size * size;
    int childOff = childIndexOff + childSize;

    float t = intersectBox(ray, origin, dims);
    if (t < 0.0) {
        return;
    }
    
    int childSubSize = readInt(array, childOff + readInt(array, childIndexOff) + 3);

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

    float childScale = scale / float(childSubSize);
    vec3 childDims = vec3(scale);

    // Step through the grid while we're inside it
    for (int i = 0; i < 3*size; i++) {
        if (c.x < 0 || c.y < 0 || c.z < 0 || c.x >= size || c.y >= size || c.z >= size) {
            return;
        }
        if (costVis) {
            ray.light.b += 0.01;
        }
        int vi = readInt(array, voxelsOff + ci);
        if (vi > 0) {
            if (costVis) {
                ray.light.b += 0.1;
            }
            int coff = childOff + readInt(array, childIndexOff + vi - 1);
            vec3 childOrigin = origin + vec3(c) * scale;
            float ct = forceIntersectBox(ray, childOrigin, childDims);
            vec3 ccp = ray.o + ray.d * (ct + childScale * 0.001);
            intersectGridLeaf(array, ray, coff, closestHit, cstep, deltaDist, childOrigin, childSubSize, childScale, ccp);
            if (closestHit.index >= 0) {
                return;
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

#else

#define MAX_HITS 16
#define MAX_TRI_HITS 8

void intersectGridLeaf(in Array array, inout Ray ray, in int headOff, inout Hit closestHit) {
    vec3 origin = readVec3(array, headOff);
    int osize = readInt(array, headOff + 3);
    int size = abs(osize);
    bool isLeaf = osize > 0;
    vec3 dims = readVec3(array, headOff + 4);
    int childSize = readInt(array, headOff + 7);
    float scale = dims.x / float(size);
    int voxelsOff = headOff + 8;
    int childIndexOff = voxelsOff + size * size * size;
    int childOff = childIndexOff + childSize;

    float t = intersectBox(ray, origin, dims);
    if (t < 0.0) {
        return;
    }

    // Map hit to voxel coordinates
    Ray tray = ray;
    tray.o = tray.o + tray.d * (t+0.001*scale);

    vec3 cf = toGrid(tray.o, scale, origin);
    ivec3 c = ivec3(cf);

    vec3 deltaDist = vec3(
        length(ray.d * ray.invD.x),
        length(ray.d * ray.invD.y),
        length(ray.d * ray.invD.z)
    );
    ivec3 cstep = ivec3(sign(ray.d));
    vec3 next = (max(vec3(0.0), vec3(cstep)) + vec3(cstep) * (vec3(c) - cf)) * deltaDist;

    int ci = c.z * size * size + c.y * size + c.x;

    // Step through the grid while we're inside it
    for (int i = 0; i < 3*size; i++) {
        if (c.x < 0 || c.y < 0 || c.z < 0 || c.x >= size || c.y >= size || c.z >= size) {
            return;
        }
        if (costVis) {
            ray.light.b += 0.01;
        }
        int vi = readInt(array, voxelsOff + ci);
        if (vi > 0) {
            if (costVis) {
                ray.light.g += 0.1;
            }
            int coff = childIndexOff + (vi - 1) * childSize;
            intersectTris(array, ray, coff, childSize, closestHit);
            ivec3 p = ivec3(toGrid(ray.o + ray.d * closestHit.distance, scale, origin));
            if (all(equal(p, c))) {
                return;
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

void intersectGridNode(in Array array, inout Ray ray, in int headOff, inout Hit closestHit) {
    vec3 origin = readVec3(array, headOff);
    int osize = readInt(array, headOff + 3);
    int size = abs(osize);
    bool isLeaf = osize > 0;
    vec3 dims = readVec3(array, headOff + 4);
    int childSize = readInt(array, headOff + 7);
    float scale = dims.x / float(size);
    int voxelsOff = headOff + 8;
    int childIndexOff = voxelsOff + size * size * size;
    int childOff = childIndexOff + childSize;

    float t = intersectBox(ray, origin, dims);
    if (t < 0.0) {
        return;
    }

    // Map hit to voxel coordinates
    Ray tray = ray;
    tray.o = tray.o + tray.d * (t+0.001*scale);

    vec3 cf = toGrid(tray.o, scale, origin);
    ivec3 c = ivec3(cf);

    vec3 deltaDist = vec3(
        length(ray.d * ray.invD.x),
        length(ray.d * ray.invD.y),
        length(ray.d * ray.invD.z)
    );
    ivec3 cstep = ivec3(sign(ray.d));
    vec3 next = (max(vec3(0.0), vec3(cstep)) + vec3(cstep) * (vec3(c) - cf)) * deltaDist;

    int ci = c.z * size * size + c.y * size + c.x;

    int hits[MAX_HITS];
    int hitCount = 0;

    // Step through the grid while we're inside it
    for (int i = 0; i < 3*size && hitCount < MAX_HITS; i++) {
        if (c.x < 0 || c.y < 0 || c.z < 0 || c.x >= size || c.y >= size || c.z >= size) {
            break;
        }
        if (costVis) {
            ray.light.b += 0.01;
        }
        int vi = readInt(array, voxelsOff + ci);
        if (vi > 0) {
            hits[hitCount++] = vi;
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

    for (int i = 0; i < hitCount; i++) {
        int vi = hits[i];
        int coff = childOff + readInt(array, childIndexOff + vi - 1);
        intersectGridLeaf(array, ray, coff, closestHit);
        if (closestHit.index >= 0) {
            return;
        }
    }

}

#endif