
// Four-level voxel grid to accelerate triangle search

uniform sampler2D arrayTex;
uniform highp isampler2D iarrayTex;
uniform highp int arrayTexWidth;

int readInt(Array array, int index) {
    int v = index / array.width;
    int u = index - (v * array.width);
    return texelFetch(iarrayTex, ivec2(u, v), 0).r;
}

float readFloat(Array array, int index) {
    int v = index / array.width;
    int u = index  - (v * array.width);
    return texelFetch(arrayTex, ivec2(u, v), 0).r;
}

vec3 readVec3(Array array, int index) {
    return vec3(
        readFloat(array, index),
        readFloat(array, index + 1),
        readFloat(array, index + 2)
    );
}

void intersectTri(in Array array, inout Ray ray, in int triIndex, inout Hit closestHit) {
    int off = 4 + triIndex * 18;
    vec3 v0 = readVec3(array, off);
    vec3 e1 = readVec3(array, off+3);
    vec3 e2 = readVec3(array, off+6);

    intersectTri(v0, e1, e2, triIndex, ray, closestHit);
}

vec3 triNormal(in Array array, in vec3 point, in int triIndex) {
    int off = 4 + triIndex * 18;

    vec3 e1 = readVec3(array, off+3);
    vec3 e2 = readVec3(array, off+6);
    
    float u = clamp(dot(point, e1), 0.0, 1.0);
    float v = clamp(dot(point, e2), 0.0, 1.0);
    vec3 n0 = readVec3(array, off+9);
    vec3 n1 = readVec3(array, off+12);
    vec3 n2 = readVec3(array, off+15);

    return normalize(mix(mix(n0, n2, v), n1, u));
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
        } else if (vi < 0) {
            if (ray.lastTested != (-vi) - 1) {
                if (costVis) {
                    ray.light.r += 0.1;
                }
                intersectTri(array, ray, (-vi) - 1, closestHit);
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

void intersectGridNodeLevel2(in Array array, inout Ray ray, in int headOff, inout Hit closestHit, 
                        ivec3 cstep, vec3 deltaDist, vec3 origin, int size, float scale, vec3 cp) {

    int childSize = readInt(array, headOff + 7);

    int voxelsOff = headOff + 8;
    int childIndexOff = voxelsOff + size * size * size;
    int childOff = childIndexOff + childSize;

    int childSubSize = abs(readInt(array, childOff + readInt(array, childIndexOff) + 3));

    vec3 cf = toGrid(cp, scale, origin);
    ivec3 c = ivec3(cf);

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

void intersectGridNodeLevel1(in Array array, inout Ray ray, in int headOff, inout Hit closestHit, 
                        ivec3 cstep, vec3 deltaDist, vec3 origin, int size, float scale, vec3 cp) {

    int childSize = readInt(array, headOff + 7);

    int voxelsOff = headOff + 8;
    int childIndexOff = voxelsOff + size * size * size;
    int childOff = childIndexOff + childSize;

    int childSubSize = abs(readInt(array, childOff + readInt(array, childIndexOff) + 3));

    vec3 cf = toGrid(cp, scale, origin);
    ivec3 c = ivec3(cf);

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
            intersectGridNodeLevel2(array, ray, coff, closestHit, cstep, deltaDist, childOrigin, childSubSize, childScale, ccp);
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

void intersectGridNodeLevel0(in Array array, inout Ray ray, in int headOff, inout Hit closestHit) {
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
    
    int childSubSize = abs(readInt(array, childOff + readInt(array, childIndexOff) + 3));

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
            intersectGridNodeLevel1(array, ray, coff, closestHit, cstep, deltaDist, childOrigin, childSubSize, childScale, ccp);
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