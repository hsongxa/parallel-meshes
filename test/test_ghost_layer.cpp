/*
  This file is part of parallel-meshes.
  parallel-meshes implements a family of unstructured 3D distributed meshes.
 
  Copyright (C) 2023  hsongxa

  This program is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#include "tests.h"
#include "parallel_mesh.h"

#include <mpi.h>
#include <vector>
#include <iostream>

int test_ghost_layer(int rank_to_display)
{
  using namespace pmh;

  int size, rank;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  if (size < 2)
  {
    std::cout << "Please run the ghost layer test with at least 2 processes!" << std::endl;
    return 1;
  }

  // hex mesh
  parallel_mesh_3d<double, hex_t> hex_mesh;
  if (rank == 0)
  {
    std::vector<double> vertices{0, 0, 0, 1, 0, 0, 2, 0, 0,
                                 0, 1, 0, 1, 1, 0, 2, 1, 0,
                                 0, 2, 0, 1, 2, 0, 2, 2, 0,
                                 0, 0, 1, 1, 0, 1, 2, 0, 1,
                                 0, 1, 1, 1, 1, 1, 2, 1, 1,
                                 0, 2, 1, 1, 2, 1, 2, 2, 1,
                                 0, 0, 2, 1, 0, 2, 2, 0, 2,
                                 0, 1, 2, 1, 1, 2, 2, 1, 2,
                                 0, 2, 2, 1, 2, 2, 2, 2, 2};
    std::vector<int> connect{0, 0, 1, 4, 3, 9, 10, 13, 12,
                             0, 1, 2, 5, 4, 10, 11, 14, 13,
                             0, 3, 4, 7, 6, 12, 13, 16, 15,
                             0, 4, 5, 8, 7, 13, 14, 17, 16,
                             0, 9, 10, 13, 12, 18, 19, 22, 21,
                             0, 10, 11, 14, 13, 19, 20, 23, 22,
                             0, 12, 13, 16, 15, 21, 22, 25, 24,
                             0, 13, 14, 17, 16, 22, 23, 26, 25};
    hex_mesh.fill_local_vertices(vertices.begin(), vertices.end());
    hex_mesh.fill_local_cells(connect.begin(), connect.end());
    hex_mesh.construct_topology();
  }

  // partition the mesh to two parts
  std::vector<int> partition_scheme;
  if (rank == 0) partition_scheme = {0, 0, 0, 0, 1, 1, 1, 1};
  hex_mesh.repartition(MPI_COMM_WORLD, partition_scheme.cbegin(), partition_scheme.cend());

  // construct ghost layer
  hex_mesh.construct_ghost_layer(MPI_COMM_WORLD);

  if (rank == rank_to_display)
  {
    std::cout << std::endl << "====== Partitioned Mesh On Rank " << rank << " ======" << std::endl;
    std::cout << "------ hex mesh vertices ------" << std::endl;
    hex_mesh.print_vertices(std::cout);
    std::cout << "------ hex mesh cell connectivities ------" << std::endl;
    hex_mesh.print_connectivity(std::cout);
    std::cout << std::endl;
  }

  return 0;
}
