// WebGL 2 shader to ray trace.

uniform sampler2D arrayTex;
uniform highp isampler2D iarrayTex;
uniform highp int arrayTexWidth;

uniform bool costVis;

uniform float deviceEpsilon;
uniform float deviceEpsilonTrace;

struct Hit {
    int index;
    float distance;
};

struct Ray {
    vec3 o;
    vec3 d;
    vec3 transmit;
    vec3 light;
    float pathLength;
    int lastTested;
};

struct Array {
    int width;
};

vec3 fromGrid(vec3 coord, float scale, vec3 origin) {
    return origin + (coord * scale);
}

vec3 toGrid(vec3 point, float scale, vec3 origin) {
    return (point - origin) / scale;
}

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

float intersectBox(in Ray ray, in vec3 origin, in vec3 dims) {
    vec3 bmin = origin;
    vec3 bmax = origin + dims;
    vec3 invD = 1.0 / ray.d;
    vec3 t0 = (bmin - ray.o) * invD;
    vec3 t1 = (bmax - ray.o) * invD;
    vec3 swaps = (sign(invD) + 1.0) * 0.5;
    vec3 swaps2 = 1.0 - swaps;
    vec3 t0_ = t0 * swaps + t1 * swaps2;
    vec3 t1_ = t0 * swaps2 + t1 * swaps;

    float tmin = max(t0_.x, max(t0_.y, t0_.z));
    float tmax = min(t1_.x, min(t1_.y, t1_.z));

    if (tmax <= tmin) {
        return -1.0;
    }
    return max(0.0, tmin);
}

void intersectTri(in Array array, inout Ray ray, in int triIndex, inout Hit closestHit) {
        int off = 4 + triIndex * 18;
        vec3 v0 = readVec3(array, off);
        vec3 e1 = readVec3(array, off+3);
        vec3 e2 = readVec3(array, off+6);

        vec3 h = cross(ray.d, e2);
        float a = dot(e1, h);

        if (a > -0.00001 && a < 0.00001) {
            return;
        }

        float f = 1.0 / a;
        vec3 s = ray.o - v0;
        float u = f * dot(s, h);

        if (u < 0.0 - deviceEpsilon || u > 1.0 + deviceEpsilon) {
            return;
        }

        vec3 q = cross(s, e1);
        float v = f * dot(ray.d, q);

        if (v < 0.0 - deviceEpsilon || u + v > 1.0 + deviceEpsilon*2.0) {
            return;
        }

        // at this stage we can compute t to find out where
        // the intersection point is on the line
        float t = f * dot(e2, q);

        if (t < 0.000001 || t > closestHit.distance - deviceEpsilon*0.1) {
            return;
        }

if (costVis) {
        if (u < 0.025 || v < 0.025 || u+v > 0.975) {
            ray.light.g += 1.0;
            ray.light.b += 1.0;
        }
}

        closestHit.index = triIndex;
        closestHit.distance = t;
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

#ifndef BREADTH_FIRST

void intersectGridLeaf(in Array array, inout Ray ray, in int headOff, inout Hit closestHit, ivec3 cstep, vec3 deltaDist) {
    vec3 origin = readVec3(array, headOff);
    int size = readInt(array, headOff + 3);
    vec3 dims = readVec3(array, headOff + 4);
    int childSize = readInt(array, headOff + 7);
    
    float scale = dims.x / float(size);

    int voxelsOff = headOff + 8;
    int childIndexOff = voxelsOff + size * size * size;
    int childOff = childIndexOff + childSize;

    float t = intersectBox(ray, origin, dims);

    vec3 cf = toGrid(ray.o + ray.d * (t + 0.001 * scale), scale, origin);
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

    vec3 cf = toGrid(ray.o + ray.d * (t + 0.001 * scale), scale, origin);
    ivec3 c = ivec3(cf);

    vec3 deltaDist = vec3(
        length(ray.d / ray.d.x),
        length(ray.d / ray.d.y),
        length(ray.d / ray.d.z)
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
                ray.light.b += 0.1;
            }
            int coff = childOff + readInt(array, childIndexOff + vi - 1);
            intersectGridLeaf(array, ray, coff, closestHit, cstep, deltaDist);
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
        length(ray.d / ray.d.x),
        length(ray.d / ray.d.y),
        length(ray.d / ray.d.z)
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
        length(ray.d / ray.d.x),
        length(ray.d / ray.d.y),
        length(ray.d / ray.d.z)
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