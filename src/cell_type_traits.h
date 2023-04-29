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

#ifndef CELL_TYPE_TRAITS_H
#define CELL_TYPE_TRAITS_H

#include "cell_operators.h"
#include "cell_topologies.h"

#include <variant>

namespace pmh {

  template<typename CT>
  struct cell_type_traits
  {
    using operator_type = single_shape_operator<CT>;
    using topology_type = single_shape_topology<CT>;
  };

  template<typename... CTS>
  struct cell_type_traits<std::variant<CTS...>>
  {
    using operator_type = mixed_shape_operator<CTS...>;
    using topology_type = mixed_shape_topology<CTS...>;
  };

} // namespace pmh

#endif
