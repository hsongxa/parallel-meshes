# parallel-meshes
A family of 3D parallel unstructured meshes implemented in modern C++ and MPI. Provided with a partitionining scheme, the implementation distributes mesh cells to MPI processes and constructs ghost layers from local cells' face neighbors.

Note that in this cimplementation, the construction of ghost layers is kept separate from the partitioning/repartitioning step. One can choose to construct the ghost layer when needed, or omit it altogether if it is not required by the numerical method. The construction of the ghost layer follows the ideas presented in

**J.M. Patchett, B. Nouanesengesy, J. Pouderoux, J. Ahrens and H. Hagen, "Parallel multi-layer ghost cell generation for distributed unstructured grids", *2017 IEEE 7th Symposium on Large Data Analysis and Visualization (LDAV)*,  Phoenix, AZ, USA, 2017, pp. 84-91**

The partitioning scheme can be provided by any partitioner, which could be graph based or space-filling curve based, such as the one from my other repository, `parallel-mesh-partitioner`, at https://github.com/hsongxa/parallel-mesh-partitioner.

The implementation includes meshes with single cell shape, such as tetrahedral mesh or hexahedral mesh, and meshes with mixed cell shapes of tetrahedron, hexahedron, pyramid, and wedge (triangular prism).

### Requirement

The code is built and tested with `g++ 9.3.0` and `MPICH 3.4.2`. For meshes with mixed cell shapes, the implementation uses the C++17 standard library's `std::variant<>` hence it should be built with at least `-std=c++17`.

### Usage

The code is all in C++ templates so there is no need to install. Simply include the header files in the `/src` directory into your projects.

Example usage of the code can be found in `/test`. To run the examples, go to the subfolder and run **`make`**. The path to the MPI installation is assumed to be `/usr/local/bin/`. For a different path, change in the **`makefile`** accordingly before running **`make`**.
