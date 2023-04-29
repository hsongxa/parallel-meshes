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

#ifndef CELL_SHAPES_H
#define CELL_SHAPES_H

namespace pmh {

  // the shape enum
  enum class shape_3d
  {
    HEXAHEDRON = 0,
    WEDGE,
    PYRAMID,
    TETRAHEDRON
  };

} // namespace pmh

#endif
