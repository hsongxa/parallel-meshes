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

#ifndef VERTEX_TYPES_H
#define VERTEX_TYPES_H

#include <tuple>
#include <utility>
#include <cassert>

namespace pmh {

  template<typename T>
  struct point_3d
  {
    T x;
    T y;
    T z;
  };

  template<typename T, typename I>
  struct vertex_3d : public point_3d<T>
  {
    I id; // vertex index local to the rank
  };

  // comparer used by kd_tree
  template<typename T, typename I>
  struct vertex_3d_comparer
  {
    int operator()(const vertex_3d<T, I>& va, const vertex_3d<T, I>& vb, int dim) const
    {
      assert(dim >= 0 && dim <= 2);
      if (dim == 0) return va.x < vb.x ? -1 : (va.x == vb.x ? 0 : 1);
      if (dim == 1) return va.y < vb.y ? -1 : (va.y == vb.y ? 0 : 1);
      if (dim == 2) return va.z < vb.z ? -1 : (va.z == vb.z ? 0 : 1);
    }
  };

  // frequently, we need to identify a cell face by its four
  // indices to vertices, hence the following utility (for
  // a triangular face, one of the indices in the tuple is
  // negative so it is always ordered to the first)

  template<typename I>
  using vid_tuple = std::tuple<I, I, I, I>;

  // order items of vid_tuple in-place
  template<typename I>
  void order_vid_tuple(vid_tuple<I>& tuple)
  {
    I& item0 = std::get<0>(tuple);
    I& item1 = std::get<1>(tuple);
    if (item0 > item1) std::swap(item0, item1);

    I& item2 = std::get<2>(tuple);
    I& item3 = std::get<3>(tuple);
    if (item2 > item3) std::swap(item2, item3);

    if (item0 > item2) std::swap(item0, item2);
    if (item1 > item3) std::swap(item1, item3);
    if (item1 > item2) std::swap(item1, item2);
  }

} // namespace pmh

#endif
