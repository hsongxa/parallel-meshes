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

#ifndef CELL_OPERATORS_H
#define CELL_OPERATORS_H

#include "cell_types.h"

#include <variant>

namespace pmh {

  // CT - the single cell type
  template<typename CT>
  struct single_shape_op
  {
    static shape_3d shape(const CT& cell) { return cell.shape; }

    static integer_type cell_vertex(const CT& cell, int v)
      { return cell._connectivity[v]; }

    static void cell_vertex(CT& cell, int v, integer_type vertex_id)
      { cell._connectivity[v] = vertex_id; }

    static hf_handle_t twin_hf(const CT& cell, int f)
      { return cell._twin_hfs[f]; }

    static void twin_hf(CT& cell, int f, hf_handle_t hf)
      { cell._twin_hfs[f] = hf; }
  };

  // CTS - multiple cell types
  template<typename... CTS>
  struct mixed_shape_op
  {
    using cell_type = std::variant<CTS...>;

    static shape_3d shape(const cell_type& cell)
      { return std::visit([](auto& c) { return c.shape; }, cell); }

    static integer_type cell_vertex(const cell_type& cell, int v)
      { return std::visit([v](auto& c) { return c._connectivity[v]; }, cell); }

    static void cell_vertex(cell_type& cell, int v, integer_type vertex_id)
      { std::visit([v, vertex_id](auto& c) { c._connectivity[v] = vertex_id; }, cell); }

    static hf_handle_t twin_hf(const cell_type& cell, int f)
      { return std::visit([f](auto& c) { return c._twin_hfs[f]; }, cell); }

    static void twin_hf(cell_type& cell, int f, hf_handle_t hf)
      { std::visit([f, hf](auto& c) { c._twin_hfs[f] = hf; }, cell); }
  };

  template<typename CT>
  struct cell_operator_traits
  {
    using operator_type = single_shape_op<CT>;
  };

  template<typename... CTS>
  struct cell_operator_traits<std::variant<CTS...>>
  {
    using operator_type = mixed_shape_op<CTS...>;
  };

} // namespace pmh

#endif
