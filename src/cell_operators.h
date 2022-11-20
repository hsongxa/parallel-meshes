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

#include <variant>

namespace pmh {

  struct single_shape_operator
  {
    template<typename CT>
    static integer_type cell_vertex(const CT& cell, int v)
      { return cell._connectivity[v]; }

    template<typename CT>
    static void cell_vertex(CT& cell, int v, integer_type vertex_id)
      { cell._connectivity[v] = vertex_id; }

    template<typename CT>
    static hf_handle_t twin_hf(const CT& cell, int f)
      { return cell._twin_hfs[f]; }

    template<typename CT>
    static void twin_hf(CT& cell, int f, hf_handle_t hf)
      { cell._twin_hfs[f] = hf; }
  };

  // CTS - multiple cell types
  struct mixed_shape_operator
  {
    template<typename... CTS>
    static integer_type cell_vertex(const std::variant<CTS...>& cell, int v)
      { return std::visit([v](auto&& c) { return c._connectivity[v]; }, cell); }

    template<typename... CTS>
    static void cell_vertex(std::variant<CTS...>& cell, int v, integer_type vertex_id)
      { std::visit([v, vertex_id](auto&& c) { c._connectivity[v] = vertex_id; }, cell); }

    template<typename... CTS>
    static hf_handle_t twin_hf(const std::variant<CTS...>& cell, int f)
      { return std::visit([f](auto&& c) { return c._twin_hfs[f]; }, cell); }

    template<typename... CTS>
    static void twin_hf(std::variant<CTS...>& cell, int f, hf_handle_t hf)
      { std::visit([f, hf](auto&& c) { c._twin_hfs[f] = hf; }, cell); }
  };

} // namespace pmh

#endif
