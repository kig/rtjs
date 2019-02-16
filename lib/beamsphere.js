/*
    The BeamCube is an acceleration structure that tesselates the search box and connects the faces, creating beams.
    Each ray intersecting with the search box is bound inside a single beam.
    If a ray enters face A and exits face B, the beam for the ray is (A, B).

    Each beam stores a list of primitives intersecting with the beam.
    The primitives are sorted by their position along the beam.
    The beam stores primitive near-far ranges for each primitive.

    To find the nearest hit primitive, you do ray-primitive intersection tests along the beam.
    Start by finding if the ray's origin intersects with any near-far range.
    If yes, test the ray against the primitives in those ranges and their intersections in the ray's direction.
    Closest hit is the closest hit in the beam and you can stop the search.
    If there's no hit, go to the next nearest range and ranges intersecting with it and repeat the test.
    Continue until you find a hit or the ray exits the beam.

    If a primitive fully covers the beam, the ray must hit it or a primitive closer to it.
    This in effect splits the beam in two, with the latter part of the beam not required for testing.
    
    In the trivial case where there are no primitives in front of the beam-covering primitive,
    you only need to intersect the beam-covering primitive.
    
    If a set of primitives fully covers the beam, you only need to intersect the primitives in the covering set.

    To accelerate covering set intersection, you create a grid of beams from the box entry face to the rear plane of the covering set.
    Each grid beam stores the primitives intersecting it.
    To intersect the covering set, first find the grid beam for the ray, then intersect with the primitives inside the beam.

    A similar strategy can be used for accelerating missing rays in case there is no covering set.
    Create beams from the box entry face to the rear plane of a partially covering set.
    Find out which beams have no intersections, then continue them to the box exit face.
    If a ray is in a no-intersection grid beam, and the beam is also no-intersection to the exit face,
    the ray doesn't intersect anything.
    If a ray is in a no-intersection grid beam, and the beam has a covering set, intersect with the covering set.
    If there is a partially covering set, continue testing against that.

    How to construct the beamcube: 
    1) Compute the bounding box of the primitives
    2) Tesselate the bounding box to faces.
    3) Create all (X,Y)-beams between the tesselated faces. For this construction, (X,Y) and (Y,X) are the same and only one is needed.
    4) Iterate through all the primitives, adding the primitive to all the beams it intersects. If the primitive covers the entire beam, add it to beam-covering sets.
    5) Sort the beam primitives by their intersection range.
    6) Bundle intersecting ranges.
    7) Find beam-covering sets of primitives in each beam.
    8) Compute coverage beamgrid: for each beam-covering set, create beamgrids at the bounding planes of the covering set.

    Computational complexity:

    Box tesselated into 6xNxN faces, each face can see 5xNxN other faces. Faces sharing their plane can't connect to each other.
    The number of beams is 5N^4 + 4N^4 + 3N^4 + 2N^4 + N^4 = 15N^4
    If N is 8, you would have 61440 beams.

    Each primitive needs to be added to all intersecting beams.
    Suppose you have a game model with 100 thousand triangles.
    You would have to do 6.144e9 triangle-beam intersections to construct the acceleration structure.
    If a primitive on average can be found in 3N^2 beams, you'd end up with 192 pointers to each primitive with N = 8.
    These would require 76.8 megabytes of storage with 4-byte pointers. 
    If you're also storing 4 bytes of range data with the primitives, the memory requirements would grow to 153 megabytes.

    The coverage acceleration beamgrid, in the case that you have two full coverage sets from each end of the beam (i.e. beam through a convex model).
    You'd have two beamgrids. If each is 32x32, that's 2048 extra beams per box search beam, each with a pointer to the intersection range.

    That'd add 61440 * 2048 * 4 bytes = 503 megabytes.

    The nice part about having a 32x32 coverage set grid is that an 8x8x32x32 BeamCube would get you avg 1 primitive intersection test per ray.

    Noodling with the idea of a 700 MB acceleration structure for a 100ktri model. 
    Would hopefully get avg 1 primitive intersection test per ray. 
    Couple with spatial path reconnection and 8 gigaray frame budget and it might do 1080p30Hz with 1000 paths per pixel.
    Add in temporal path reconnection, denoising, importance sampling based on path cache energies, and you should be getting very high quality interactives.

    vec2 signNotZero(vec2 v) {
        return vec2((v.x >= 0.0) ? +1.0 : -1.0, (v.y >= 0.0) ? +1.0 : -1.0);
    }

    vec2 toOct(in vec3 v) {
        // Project the sphere onto the octahedron, and then onto the xy plane
        // &
        // Reflect the folds of the lower hemisphere over the diagonals

        vec2 p = v.xy * (1.0 / (abs(v.x) + abs(v.y) + abs(v.z)));
        return (v.z <= 0.0) ? ((1.0 - abs(p.yx)) * signNotZero(p)) : p;
    }
    
    vec3 fromOct(vec2 e) {
        vec3 v = vec3(e.xy, 1.0 - abs(e.x) - abs(e.y));
        if (v.z < 0) v.xy = (1.0 - abs(v.yx)) * signNotZero(v.xy);
        return normalize(v);
    }
*/
class Beam {
    static cmp(a, b) {
        return a.closestDistance - b.closestDistance;
    }

    constructor(startX, startY, endX, endY, sphereOrigin, sphereRadius, radius) {
        this.startX = startX;
        this.startY = startY;
        this.endX = endX;
        this.endY = endY;
        this.sphereRadius = sphereRadius;
        this.radius = radius;
        this.diskSize = 4;
        this.triangles = [];
        this.lookup = [];
        for (let i=0; i < this.diskSize*this.diskSize; i++) {
            this.lookup.push([]);
        }
        this.startPoint = add(sphereOrigin, mulS(fromOct(this.startX, this.startY), this.sphereRadius));
        this.endPoint = add(sphereOrigin, mulS(fromOct(this.endX, this.endY), this.sphereRadius));
        this.ray = {o: this.startPoint, d: normalize(sub(this.endPoint, this.startPoint))};
        this.planeUVW = orthoBasis(this.ray.d);
    }

    addTriangle(tri) {
        const startPoint = this.startPoint;
        const dir = this.ray.d;
        const d0 = sub(tri._vertices[0], startPoint);
        const d1 = sub(tri._vertices[1], startPoint);
        const d2 = sub(tri._vertices[2], startPoint);
        const closestDistance = min(min(dot(d0, dir), dot(d1, dir)), dot(d2, dir));
        const triObj = {obj: tri, closestDistance: closestDistance};
        this.triangles.push(triObj);
    }

    sort() {
        this.triangles.sort(Beam.cmp);
    }

    constructBeamGrid() {
        if (this.triangles.length < 5) {
            return;
        }
        const uvw = this.planeUVW;
        const n = this.diskSize;
        const n1 = n - 1;
        const radius = this.radius;
        const bradius = radius * (n1/n) * 0.5;
        const tris = this.triangles;
        for (let y = 0; y < n; y++) {
            for (let x = 0; x < n; x++) {
                // const uf = ((x/n1) - 0.5) * 2;
                // const vf = ((y/n1) - 0.5) * 2;
                // const o = diskUVToPoint(this.startPoint, uvw, uf, vf, bradius);
                // const r = new Ray(o, this.ray.d);
                // const lookupTris = [];
                // for (let i = 0; i < tris.length; i++) {
                //     if (tris[i].obj.intersectCylinder(r, bradius/n1)) {
                //         // this.lookup[y*n + x] = tris[i];
                //         lookupTris.push(tris[i]);
                //         break;
                //     }
                // }
                // // console.log(lookupTris.length, tris.length);
                // this.lookup[y*n + x] = lookupTris;

                // create beam N^2 grid
                const uf0 = ((x/n1) - 0.5) * 2;
                const vf0 = ((y/n1) - 0.5) * 2;
                const o = diskUVToPoint(this.startPoint, uvw, uf0, vf0, bradius);
                const r = new Ray(o, vec3(0));
                const lookup = this.lookup[y*n + x];
                for (let ty = 0; ty < n; ty++) {
                    for (let tx = 0; tx < n; tx++) {
                        lookup[ty*n + tx] = null;
                        if (x !== tx || y !== ty) {
                            const uf1 = ((tx/n1) - 0.5) * 2;
                            const vf1 = ((ty/n1) - 0.5) * 2;
                            const p = diskUVToPoint(this.endPoint, uvw, uf1, vf1, bradius);
                            r.d = normalize(sub(p, o));
                            // const lookupTris = [];
                            tris: 
                            for (let i = 0; i < tris.length; i++) {
                                if (tris[i].obj.intersect(r)) {
                                    lookup[ty*n + tx] = tris[i].obj;
                                    break tris;
                                }
                            }
                            // lookup[ty*n + tx] = lookupTris;
                            // this.allBeams.push(lookupTris);
                        }
                    }
                }
            }
        }
    }

    fastIntersect(ray, rayStart, rayEnd) {
        if (this.triangles.length < 5) {
            return this.intersect(ray);
        }
        const entryPoint = add(ray.o, mulS(ray.d, rayStart));
        const exitPoint = add(ray.o, mulS(ray.d, rayEnd));
        const n = this.diskSize;
        const n1 = n-1;
        const bradius = this.radius * (n1/n) * 0.5;
        const uv0 = diskPointToUV(entryPoint, this.startPoint, this.planeUVW, bradius);
        const x0 = Math.floor((uv0.u / 2 + 0.5) * n1 - 1e-9);
        const y0 = Math.floor((uv0.v / 2 + 0.5) * n1 - 1e-9);
        const uv1 = diskPointToUV(exitPoint, this.endPoint, this.planeUVW, bradius);
        const x1 = Math.floor((uv1.u / 2 + 0.5) * n1 - 1e-9);
        const y1 = Math.floor((uv1.v / 2 + 0.5) * n1 - 1e-9);
        const off0 = y0*n + x0;
        const off1 = y1*n + x1;
        // if (this.lookup[off0][off1]) {
        //     const lookup = this.lookup[off0][off1];
        //     const closestHit = {obj: null, distance: Infinity};
        //     for (let i = 0; i < lookup.length; i++) {
        //         const tri = lookup[i];
        //         if (tri.closestDistance >= closestHit.distance) {
        //             return closestHit;
        //         }
        //         BeamSphere.primitiveTests++;
        //         const hit = tri.obj.intersect(ray);
        //         if (hit && hit.distance < closestHit.distance) {
        //             closestHit.obj = hit.obj;
        //             closestHit.distance = hit.distance;
        //         }
        //     }
        //     return closestHit.obj ? closestHit : null;
        // }
        // return null;

        // let x = x0;
        // let y = y0;
        // let dx = x1-x0;
        // let dy = y1-y0;
        // let de = abs(dy / dx);
        // let e = 0;
        // let sdy = dy >= 0 ? 1 : -1;
        // let sdx = dx >= 0 ? 1 : -1;
        // while (y !== y1 && x !== x1) {
        //     const off = y*n + x;
        //     if (this.lookup[off]) {
        //         const tris = this.lookup[off];
        //         for (let i = 0; i < tris.length; i++) {
        //             BeamSphere.primitiveTests++;
        //             const t = tris[i].obj.intersect(ray);
        //             if (t) {
        //                 return t;
        //             }
        //         }
        //     }
        //     e += de;
        //     if (e >= 0.5) {
        //         y = y + sdy;
        //         e -= 1;
        //     } else {
        //         x = x + sdx;
        //     }
        // }
        if (this.lookup[off0] && this.lookup[off0][off1]) {
            BeamSphere.primitiveTests++;
            return this.lookup[off0][off1].intersectAsPlane(ray);
        }
        // if (this.lookup[off0] && this.lookup[off0][0]) {
        //     BeamSphere.primitiveTests++;
        //     return this.lookup[off0][0].obj.intersectAsPlane(ray);
        // }
        return null;
    }

    intersectTriangle(tri) {
        return tri.intersectCylinder(this.ray, this.radius);
    }

    intersect(ray) {
        const tris = this.triangles;
        let hit = {obj: null, distance: Infinity};
        for (let i = 0; i < tris.length; i++) {
            const tri = tris[i];
            if (tri.closestDistance >= hit.distance) {
                break;
            }
            const triHit = tri.obj.intersect(ray);
            BeamSphere.primitiveTests++;
            if (triHit && triHit.distance < hit.distance) {
                hit = triHit;
            }
        }
        return hit.obj ? hit : null;
    }
}

class BeamSphere {

    constructor(origin, radius, beamsPerSide) {
        this.origin = origin;
        this.sphere = new Sphere(origin, radius);
        this.beamsPerSide = beamsPerSide;
        const beams = this.beams = [];
        this.allBeams = [];
        for (let y = 0; y < beamsPerSide; y++) {
            for (let x = 0; x < beamsPerSide; x++) {
                beams[y*beamsPerSide + x] = [];
            }
        }
        const d = 1/beamsPerSide * 0.5;
        for (let y = 0; y < beamsPerSide; y++) {
            for (let x = 0; x < beamsPerSide; x++) {
                for (let ty = 0; ty < beamsPerSide; ty++) {
                    for (let tx = 0; tx < beamsPerSide; tx++) {
                        if (x !== tx || y !== ty) {
                            const beam = new Beam(
                                x/beamsPerSide+d, y/beamsPerSide+d, 
                                tx/beamsPerSide+d, ty/beamsPerSide+d, 
                                origin, radius, 
                                1.2*Math.PI*radius/beamsPerSide
                            );
                            beams[y*beamsPerSide + x].push(beam);
                            this.allBeams.push(beam);
                        }
                    }
                }
            }
        }
    }

    sort() {
        this.allBeams.forEach(b => {
            b.sort();
            b.constructBeamGrid();
        });
    }

    add(obj) {
        for (let i = 0; i < this.allBeams.length; i++) {
            if (this.allBeams[i].intersectTriangle(obj)) {
                this.allBeams[i].addTriangle(obj);
                BeamSphere.addedTriangles++;
            }
        }
    }

    intersect(ray) {
        const t = this.sphere.intersectBoth(ray);
        BeamSphere.sphereTests++;
        if (!t) {
            return t;
        }

        const beamsPerSide = this.beamsPerSide;

        const origin = this.sphere.center;

        let entry = toOct(sub(add(ray.o, mulS(ray.d, t.entryDistance)), origin));
        let exit = toOct(sub(add(ray.o, mulS(ray.d, t.exitDistance)), origin));

        let x = floor(entry.x * beamsPerSide - 1e-9);
        let y = floor(entry.y * beamsPerSide - 1e-9);
        let tx = floor(exit.x * beamsPerSide - 1e-9);
        let ty = floor(exit.y * beamsPerSide - 1e-9);

        const beam = this.beams[y * beamsPerSide + x][ty * beamsPerSide + tx];
        BeamSphere.beamTests++;
        if (!beam) {
            return null;
        }
        if (false) {
            ray.light.x = beam.startX;
            ray.light.y = beam.startY + 5 * beam.triangles.length === 0;
            ray.light.z = beam.endY + beam.endX;
            const sd = length(sub(add(ray.o, mulS(ray.d, t.entryDistance)), beam.startPoint));
            const ed = length(sub(add(ray.o, mulS(ray.d, t.exitDistance)), beam.endPoint));
            if (sd < 0.1) {
                ray.light = add(ray.light, vec3(0,1,0));
            } else if (2*(sd / beam.radius) % 1 < 0.1) {
                ray.light = add(ray.light, vec3(0,1,1));
            } else if (rayLineDistance(ray, beam.startPoint, beam.endPoint) < beam.radius*0.05) {
                ray.light = add(ray.light, vec3(0,0,1));
            } else if (ed < 0.1) {
                ray.light = add(ray.light, vec3(1,0,0));
            }
            ray.finished = true;
            return;
        }
        // return beam.intersect(ray);
        return beam.fastIntersect(ray, t.entryDistance, t.exitDistance);
    }

}

BeamSphere.sphereTests = 0;
BeamSphere.beamTests = 0;
BeamSphere.primitiveTests = 0;
BeamSphere.addedTriangles = 0;