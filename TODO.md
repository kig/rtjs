### TODO

[] Lookdev an artwork look in Photoshoppe
    [] Diamonds, stripes
    [] Animations, interactivity
    [] Text, text model

[] 8x8 dense layout
    [] Create BitArray(8^3 + 8^3 * 8^3) for the grid walk
        [] Top-level VG is the first 8^3 bits
        [] VG[x,y,z] is located at bitArray[8^3 + (z * 8 * 8 + y * 8 + x) * 8^3]
    [] Create UInt16Array(8^3 * 8^3 * 16) for voxels
    [] Add triangle indices directly to voxel cells (first zero cell)
    [] Encode triangles as vertex + oct-encoded edges * radius + oct-encoded normal
    => 32 kB bit array, 8 MB voxel array, 128 kB tri array @ 32 bytes per tri
    - 4 L2 accesses (tris), 2 RAM accesses (voxels) => 4 * 125ns + 2 * 225ns = 950ns

    Sparse layout is 192 VGs, 13887 voxels, 44k prims in voxels
    => 12 kB bit array, 400 kB subgrid pointers, 176kB tri indices, 128 kB tri array @ 32 B/tri
    => Fits in L2...
    - 4 L2 for tris, 3 L2 for VG, 2x2 L2 for voxels => 11 * 125ns = 1350ns

    Keep first level dense, subgrids sparse
    - 4 L2 for tris, 2x2 L2 for voxels = 8 * 125 = 1000ns

    Anyway it's at least 4 L2s since the tris don't fit in L1, so the floor is 500ns.
    Probably best to stuff voxels into L2 so that fetching them doesn't push the tris out.
    

[] Load models with THREE ObjLoader (Sponza, etc.)
[] Deal with models with textures and UVs
[] Crazy diamonds
[] Fix bbox-tri problems in dragon.obj

[] Materials
    [] Transparent materials
    [] Plastic materials
    [] Subsurface scattering
    [] Coated materials

[] Optimize VoxelGrid representation
    [] Bytes for finding matching nodes
    [] Bitfield for finding matching nodes
    [] Grid shape defined in header

[] Multiple objects
    [] Instancing

[] Animated objects
    [] Rigid animation (transform acceleration structure)

[] Better paths
    [] BPT
    [] VCM
    [] Path reconnection

[] Denoise
    [] Spatial denoise
    [] Temporal denoise
