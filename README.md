# parallel-meshes
A family of 3D parallel unstructured meshes implemented in modern C++ and MPI. The parallelism, i.e., partitioning, re-partitioning, and load balancing, is based on space filling curves as described in the paper

**C. Burstedde and J. Holke, "Coarse mesh partitioning for tree-based AMR", *SIAM J. Sci. Comput.*, 39(5), C364-C392, 2017**

The implementation includes meshes with single cell shape, such as tetrahedral mesh and hexahedral mesh, and meshes with mixed cell shapes of tetrahedron, hexahedron, pyramid, and wedge (triangular prism). Mesh adaptation, i.e., refinement and coarsening, is also supported.

### Requirement

The code is built and tested with `g++ 9.3.0` and `MPICH 3.4.2`.

### Usage

Example usage of the code can be found in `/examples`. To run the examples, go to the subfolder and run **`make`**. The path to the MPI installation is assumed to be `/usr/local/bin/`. For a different path, change in the **`makefile`** accordingly before running **`make`**.
