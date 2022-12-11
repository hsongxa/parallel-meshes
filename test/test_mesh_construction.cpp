/*
  Copyright (C) 2022 Hao Song

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

#include <tuple>
#include <vector>
#include <variant>
#include <iostream>

int test_mesh_construction()
{
  using namespace pmh;

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
  // hex mesh
  parallel_mesh_3d<double, hex_t> hex_mesh(0);
  hex_mesh.fill_local_vertices(vertices.begin(), vertices.end());
  hex_mesh.fill_local_cells(connect.begin(), connect.end());
  hex_mesh.construct_topology();

  std::cout << "------ hex mesh cell connectivities ------" << std::endl;
  hex_mesh.print_connectivity(std::cout);
  std::cout << "------ hex mesh twin half-facets ------" << std::endl;
  hex_mesh.print_twin_hfs(std::cout);
  std::cout << std::endl;

  std::tuple<double, double, double, double, double, double> bbox = hex_mesh.local_bounding_box();
  std::cout << "local bounding box: [" << std::get<0>(bbox) << ", " << std::get<1>(bbox) << ", ";
  std::cout << std::get<2>(bbox) << ", " << std::get<3>(bbox) << ", " << std::get<4>(bbox) << ", ";
  std::cout << std::get<5>(bbox) << "]" << std::endl;

  std::cout << "------ hex mesh in Gmsh format ------" << std::endl;
  hex_mesh.export_to_msh_format(std::cout);

  // add one more vertex and a pyramid to the list
  vertices.push_back(1);
  vertices.push_back(1);
  vertices.push_back(3);
  connect.push_back(2);
  connect.push_back(18);
  connect.push_back(19);
  connect.push_back(22);
  connect.push_back(21);
  connect.push_back(27);

  // mixed shape mesh of hex, tet, wedge, and pyramid
  parallel_mesh_3d<double, std::variant<hex_t, tet_t, wdg_t, prm_t>> mixed_mesh(0);
  mixed_mesh.fill_local_vertices(vertices.begin(), vertices.end());
  mixed_mesh.fill_local_cells(connect.begin(), connect.end());
  mixed_mesh.construct_topology();

  std::cout << std::endl;
  std::cout << "------ mixed mesh cell connectivities ------" << std::endl;
  mixed_mesh.print_connectivity(std::cout);
  std::cout << "------ mixed mesh twin half-facets ------" << std::endl;
  mixed_mesh.print_twin_hfs(std::cout);
  std::cout << std::endl;

  bbox = mixed_mesh.local_bounding_box();
  std::cout << "local bounding box: [" << std::get<0>(bbox) << ", " << std::get<1>(bbox) << ", ";
  std::cout << std::get<2>(bbox) << ", " << std::get<3>(bbox) << ", " << std::get<4>(bbox) << ", ";
  std::cout << std::get<5>(bbox) << "]" << std::endl;

  // index of local cell will be automatically converted
  // to cell handle - this applies to local cell only
  std::tuple<double, double, double> centroid = mixed_mesh.cell_centroid(8);
  std::cout << "centroid of the last cell: (" << std::get<0>(centroid);
  std::cout << ", " << std::get<1>(centroid) << ", " << std::get<2>(centroid) << ")" << std::endl;

  return 0;
}
