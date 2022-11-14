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

#include <variant>

int test_mesh_construction()
{
  using namespace pmh;

  // tet mesh
  parallel_mesh_3d<float, tet_t> tet_mesh(0);

  // hex mesh
  parallel_mesh_3d<double, hex_t> hex_mesh(0);

  // mixed shape mesh of tet, hex, and wedge
  parallel_mesh_3d<float, std::variant<tet_t, hex_t, wdg_t>> mixed_3_mesh(0);

  // mixed shape mesh of hex, tet, wedge, and pyramid
  parallel_mesh_3d<double, std::variant<hex_t, tet_t, wdg_t, prm_t>> mixed_4_mesh(0);

  return 0;
}
