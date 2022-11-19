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

#ifndef CELL_TYPES_H
#define CELL_TYPES_H

#include "cell_topologies.h"

#include <cstdint>

namespace pmh {

  using integer_type = std::int64_t;

  // the first four bits are reserved to encode the refinement level (maximum 15), and
  // the last bit (i.e., the least significant bit) is used to distinguish local vs ghost
  // and the remaining bits define the index inside a MPI rank
  //
  // half-facet is represented by a cell handle and a local face index using three bits;
  // so practically the upper limit for the number of cells is
  //
  //   2^(64 - 4 - 3 - 1 - 1) = 3.6e16
  //               ^       ^
  //             hf bits  sign bit
  //
  // per level per rank when using 64-bit integer
  //
  // value of -1 represents an invalid cell handle
  struct cell_handle_t
  {
    cell_handle_t() : _int_rep(-1) {}

    cell_handle_t(integer_type index, bool is_ghost = false, int level = 0)
      : _int_rep((index << 1) | ((integer_type)level << 56)) { if (is_ghost) ++_int_rep; }

    bool is_valid() const { return _int_rep >= 0; }

    int level() const { return (_int_rep & (integer_type)1080863910568919040) >> 56; }

    integer_type local_index() const { return (_int_rep >> 1) & (integer_type)36028797018963967; }

    bool is_ghost() const { return _int_rep % 2; }

    // can accept negative offset which means backward
    cell_handle_t forward(int offset) const
      { return cell_handle_t(local_index() + offset, level(), is_ghost()); };

    // make the underlying integer representation public
    // so it can be used in message passing
    integer_type _int_rep;
  };

  inline bool operator<(const cell_handle_t& x, const cell_handle_t& y) { return x._int_rep < y._int_rep; }
  inline bool operator==(const cell_handle_t& x, const cell_handle_t& y) { return x._int_rep == y._int_rep; }


  // half-facet is encoded by the cell's handle and the local face
  // index (using the three least significant bits)
  //
  // default to -1 to indicate invalid half-facet
  struct hf_handle_t
  {
    hf_handle_t() : _int_rep(-1) {}

    hf_handle_t(cell_handle_t cid, int face) : _int_rep(cid._int_rep << 3 | face) {}

    hf_handle_t(integer_type cid, bool is_ghost, int level, int face)
      : hf_handle_t(cell_handle_t(cid, is_ghost, level), face) {}

    bool is_valid() const { return _int_rep >= 0; }

    cell_handle_t cell_handle() const { cell_handle_t tmp; tmp._int_rep = _int_rep >> 3; return tmp; }

    int face_index() const { return _int_rep & 7; }

    // make the underlying integer representation public
    // so it can be used in message passing
    integer_type _int_rep;
  };

  // hexahedron
  struct hex_t
  {
    hex_t() : _connectivity{ -1, -1, -1, -1, -1, -1, -1, -1 } {}

    integer_type _connectivity[8];
    hf_handle_t  _twin_hfs[6];

    using topology = hex_topology;
  };

  // wedge
  struct wdg_t
  {
    wdg_t() : _connectivity{ -1, -1, -1, -1, -1, -1 } {}

    integer_type _connectivity[6];
    hf_handle_t  _twin_hfs[5];

    using topology = wdg_topology;
  };

  // pyramid
  struct prm_t
  {
    prm_t() : _connectivity{ -1, -1, -1, -1, -1 } {}

    integer_type _connectivity[5];
    hf_handle_t  _twin_hfs[5];

    using topology = prm_topology;
  }; 

  // tetrahedron
  struct tet_t
  {
    tet_t() : _connectivity{ -1, -1, -1, -1 } {}

    integer_type _connectivity[4];
    hf_handle_t  _twin_hfs[4];

    using topology = tet_topology;
  };

} // namespace pmh

#endif
