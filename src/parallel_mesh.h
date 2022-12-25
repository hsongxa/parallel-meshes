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

#ifndef PARALLEL_MESH_H
#define PARALLEL_MESH_H

#include "vertex_types.h"
#include "cell_types.h"
#include "cell_type_traits.h"
#include "mpi_datatype_traits.h"
#include "kd_tree.h"

#include <mpi.h>
#include <cstddef>
#include <limits>
#include <tuple>
#include <vector>
#include <set>
#include <map>
#include <cassert>
#include <algorithm>
#include <numeric>
#include <iostream>

namespace pmh {

  // Unstructured 3D distributed meshes with any combination
  // of cell types of tet, hex, wedge, or pyramid.
  //
  // R - floating point type for coordinates, CT - cell type
  // which could be single cell type or variant of multiple
  // cell types, OP - cell operator type that operates on CT,
  // TP - cell topology corresponding to CT
  //
  // Passing a tet_t type to CT will make a tetrahedral mesh,
  // a hex_t for hexahedral mesh, ..., etc. For mixed
  // tet-wdg-hex mesh, CT is std::variant<hex_t, wdg_t, tet_t>.
  // Similarly, std::variant<hex_t, wdg_t, prm_t, tet_t> will
  // give a mesh with mixed cell shapes of all four.
  //
  // The partitioning/repartitioning and the construction of
  // ghost layers are kept separate, i.e., partitioning only
  // (re)distributes local cells, not ghost cells. And the
  // construction of ghost layers should be performed after
  // the partitioning/repartitioning.
  //
  // Vertices and cells are identified by their local indices
  // on the rank - there are no global ID's for them.
  template< class R, class CT,
            class OP = typename cell_type_traits<CT>::operator_type,
            class TP = typename cell_type_traits<CT>::topology_type >
  class parallel_mesh_3d
  {
  public:
    using point_type = point_3d<R>;
    using cell_type = CT;
    using cell_op = OP;
    using cell_tp = TP;
    using size_type = std::size_t;

    // ------ local construction operations ------

    // population of vertices used by local cells (no ghosts)
    template<class InputIt>
    void fill_local_vertices(InputIt first, InputIt last);

    // population of local cell connectivities
    template<class InputIt>
    void fill_local_cells(InputIt first, InputIt last);

    // ------ MPI operations ------

    // partitioning/repartitioning
    // TODO: expand to handle cell-based property re-distribution
    template<class InputIt>
    void repartition(MPI_Comm comm, InputIt first, InputIt last);

    // one layer of ghost cells that are face neighbors of local cells
    // TODO: expand to also distribute cell-based properties for the ghost layer
    void construct_ghost_layer();

    // ------ queries of vertices, faces, and cells ------
    // ------ bounding box and centroid are needed by SFC based mesh partitioners ------

    size_type num_vertices() const { return _vertices.size(); }

    size_type num_local_cells() const { return _local_cells.size(); }

    size_type num_ghost_cells() const { return _ghost_cells.size(); }

    const point_type& get_vertex(integer_type vertex_id) const { return _vertices[vertex_id]; }

    point_type& get_vertex(integer_type vertex_id) { return _vertices[vertex_id]; }

    bool is_boundary_face(cell_handle_t cell_handle, int f) const
    {
      hf_handle_t hf = cell_op::twin_hf(get_cell(cell_handle), f);
      return !hf.is_valid();
    }

    const cell_type& get_cell(cell_handle_t cell_handle) const
    {
      assert(cell_handle.is_valid());
      if (cell_handle.is_ghost()) return _ghost_cells[cell_handle.local_index()];
      return _local_cells[cell_handle.local_index()];
    }

    cell_type& get_cell(cell_handle_t cell_handle)
    {
      assert(cell_handle.is_valid());
      if (cell_handle.is_ghost()) return _ghost_cells[cell_handle.local_index()];
      return _local_cells[cell_handle.local_index()];
    }

    std::tuple<R, R, R, R, R, R> local_bounding_box() const;

    std::tuple<R, R, R> cell_centroid(cell_handle_t cell_handle) const;

    // ------ queries of topological relationships ------

    // from facet-to-facet and vertex-to-facet relationships, all
    // other topological relationships can be derived for the given
    // cell types/topologies
    void construct_topology() { fill_twin_hfs(); fill_vtx_to_hfs(); }

    // invalid half-facet handle is returned for boundary faces
    hf_handle_t get_face_neighbor(cell_handle_t cell_handle, int f) const
      { return cell_op::twin_hf(get_cell(cell_handle), f); }

    // return any facet incident to the given vertex, priority given to boundary facet
    hf_handle_t get_incident_face(integer_type vertex_id) const
      { return _vtx_to_hfs[vertex_id]; }

    static shape_3d cell_shape(const cell_type& cell) { return cell_tp::shape(cell); }

    // TODO: other topological relationships...

    // ------ static getters for a given cell ------

    static integer_type cell_vertex(const cell_type& cell, int v)
      { return cell_op::cell_vertex(cell, v); }

    static hf_handle_t twin_hf(const cell_type& cell, int f)
      { return cell_op::twin_hf(cell, f); }

  private:
    // ------ static setters for a given cell ------

    static void cell_vertex(cell_type& cell, int v, integer_type vertex_id)
      { cell_op::cell_vertex(cell, v, vertex_id); }

    static void twin_hf(cell_type& cell, int f, hf_handle_t hf)
      { cell_op::twin_hf(cell, f, hf); }

    // ------ topology construction (involving ghost cells) ------

    void fill_twin_hfs();

    void fill_vtx_to_hfs();

    // remove duplicate vertices and update connectivities
    void remove_duplicate_vertices(std::vector<point_type>& vertices, std::vector<cell_type>& cells);

    // ------ MPI communications ------

    // it0: beginning of send scheme; it1: beginning of indices of cells to send;
    // vtx_it: output of received vertices; cell_it: output of received cell connectivities
    template<typename RandIt, typename OutIt1, typename OutIt2>
    void send_recv_cells(MPI_Comm comm, RandIt it0, RandIt it1, OutIt1 vtx_it, OutIt2 cell_it);

  private:
    std::vector<point_type>  _vertices;
    std::vector<hf_handle_t> _vtx_to_hfs;

    std::vector<cell_type>   _local_cells;
    std::vector<cell_type>   _ghost_cells;

    static constexpr int MAX_NUM_VERTICES_PER_FACE = 4;
    static constexpr int MAX_NUM_FACES_INCIDENT_TO_VERTEX = 4; // pyramid has 4, all other types have 3

  public: // mainly for debugging
    void print_connectivity(std::ostream& out) const;
    void print_twin_hfs(std::ostream& out) const;
    // export: local cells only (config=0); ghost cells only (config=1); or both local and ghost cells (config=2)
    void export_to_msh_format(std::ostream& out, int config = 0) const;
  };

  template<class R, class CT, class OP, class TP> template<class InputIt>
  void parallel_mesh_3d<R, CT, OP, TP>::fill_local_vertices(InputIt first, InputIt last)
  {
    _vtx_to_hfs.clear();
    _vertices.clear();

    while (first != last)
      _vertices.emplace_back(point_type{*first++, *first++, *first++});
  }

  template<class R, class CT, class OP, class TP> template<class InputIt>
  void parallel_mesh_3d<R, CT, OP, TP>::fill_local_cells(InputIt first, InputIt last)
  {
    _ghost_cells.clear();
    _local_cells.clear();

    while (first != last) // format of each cell is: <shape> <v1, v2, ...>
    {
      shape_3d shape = static_cast<shape_3d>(*first++);
      _local_cells.push_back(cell_op::create_cell(shape));
      cell_type& cell = _local_cells[_local_cells.size() - 1];
      for (int i = 0; i < cell_tp::num_vertices(cell); ++i)
        cell_vertex(cell, i, *first++);
    }
  }

  template<class R, class CT, class OP, class TP>
  void parallel_mesh_3d<R, CT, OP, TP>::fill_twin_hfs()
  {
    int face_vertices[MAX_NUM_VERTICES_PER_FACE];
    std::map<vid_tuple<integer_type>, hf_handle_t> dict;

    for (size_type c = 0; c < _local_cells.size(); ++c)
    {
      cell_type& cell = _local_cells[c];
      cell_handle_t cell_handle = cell_handle_t(c, false);
      for (int f = 0; f < cell_tp::num_faces(cell); ++f)
      {
        cell_tp::vertices_of_face(cell, f, face_vertices);
        int num_vertices = cell_tp::num_vertices_of_face(cell, f);

        vid_tuple<integer_type> vtuple = num_vertices < 4 ?
          std::make_tuple(integer_type(-1), cell_vertex(cell, face_vertices[0]),
                          cell_vertex(cell, face_vertices[1]), cell_vertex(cell, face_vertices[2])) :
          std::make_tuple(cell_vertex(cell, face_vertices[0]), cell_vertex(cell, face_vertices[1]),
                          cell_vertex(cell, face_vertices[2]), cell_vertex(cell, face_vertices[3]));
        order_vid_tuple(vtuple);

        auto it_twin = dict.find(vtuple);
        if (it_twin != dict.end())
        {
          hf_handle_t hf = it_twin->second;
          twin_hf(cell, f, hf);
          twin_hf(get_cell(hf.cell_handle()), hf.face_index(), hf_handle_t(cell_handle, f));
          dict.erase(it_twin);
        }
        else
        {
          twin_hf(cell, f, hf_handle_t()); // boundary facet indicated by invalid twin hf
          dict.insert(std::make_pair(vtuple, hf_handle_t(cell_handle, f)));
        }
      }
    }

    for (size_type c = 0; c < _ghost_cells.size(); ++c)
    {
      cell_type& cell = _ghost_cells[c];
      cell_handle_t cell_handle = cell_handle_t(c, true);
      for (int f = 0; f < cell_tp::num_faces(cell); ++f)
      {
        cell_tp::vertices_of_face(cell, f, face_vertices);
        int num_vertices = cell_tp::num_vertices_of_face(cell, f);

        vid_tuple<integer_type> vtuple = num_vertices < 4 ?
          std::make_tuple(integer_type(-1), cell_vertex(cell, face_vertices[0]),
                          cell_vertex(cell, face_vertices[1]), cell_vertex(cell, face_vertices[2])) :
          std::make_tuple(cell_vertex(cell, face_vertices[0]), cell_vertex(cell, face_vertices[1]),
                          cell_vertex(cell, face_vertices[2]), cell_vertex(cell, face_vertices[3]));
        order_vid_tuple(vtuple);

        auto it_twin = dict.find(vtuple);
        if (it_twin != dict.end())
        {
          hf_handle_t hf = it_twin->second;
          twin_hf(cell, f, hf);
          twin_hf(get_cell(hf.cell_handle()), hf.face_index(), hf_handle_t(cell_handle, f));
          dict.erase(it_twin);
        }
        else
        {
          twin_hf(cell, f, hf_handle_t()); // boundary facet indicated by invalid twin hf
          dict.insert(std::make_pair(vtuple, hf_handle_t(cell_handle, f)));
        }
      }
    }
  }

  // ASSUMPTION: this is always called after fill_twin_hfs()
  template<class R, class CT, class OP, class TP>
  void parallel_mesh_3d<R, CT, OP, TP>::fill_vtx_to_hfs()
  {
    int incident_faces[MAX_NUM_FACES_INCIDENT_TO_VERTEX];
    int face_vertices[MAX_NUM_VERTICES_PER_FACE];

    // allocate and initialize with invalid hfs
    _vtx_to_hfs.resize(_vertices.size());
    std::fill(_vtx_to_hfs.begin(), _vtx_to_hfs.end(), hf_handle_t());

    // first, go through ghost cells
    for (size_type c = 0; c < _ghost_cells.size(); ++c)
    {
      cell_type& cell = _ghost_cells[c];
      cell_handle_t cell_handle = cell_handle_t(c, true);

      // populate all its vertices with the first incident face
      for (int v = 0; v < cell_tp::num_vertices(cell); ++v)
      {
        integer_type vid = cell_vertex(cell, v);
        if (!_vtx_to_hfs[vid].is_valid())
        {
          cell_tp::incident_faces(cell, v, incident_faces);
          _vtx_to_hfs[vid] = hf_handle_t(cell_handle, incident_faces[0]);
        }
      }
      
      // then give priority to boundary faces
      for (int f = 0; f < cell_tp::num_faces(cell); ++f)
        if (is_boundary_face(cell_handle, f))
        {
          cell_tp::vertices_of_face(cell, f, face_vertices);
          for (int v = 0; v < cell_tp::num_vertices_of_face(cell, f); ++v)
            _vtx_to_hfs[cell_vertex(cell, face_vertices[v])] = hf_handle_t(cell_handle, f);
        }
    }

    // then, local cells - this may overwrite the info filled above,
    // which means local cell' faces get higher priority
    for (size_type c = 0; c < _local_cells.size(); ++c)
    {
      cell_type& cell = _local_cells[c];
      cell_handle_t cell_handle = cell_handle_t(c, false);

      // populate all its vertices with the first incident face
      for (int v = 0; v < cell_tp::num_vertices(cell); ++v)
      {
        integer_type vid = cell_vertex(cell, v);
        hf_handle_t hf = _vtx_to_hfs[vid];
        if (!hf.is_valid() || hf.cell_handle().is_ghost()) // overwrite ghost's hf
        {
          cell_tp::incident_faces(cell, v, incident_faces);
          _vtx_to_hfs[vid] = hf_handle_t(cell_handle, incident_faces[0]);
        }
      }
      
      // then give priority to boundary faces
      for (int f = 0; f < cell_tp::num_faces(cell); ++f)
        if (is_boundary_face(cell_handle, f))
        {
          cell_tp::vertices_of_face(cell, f, face_vertices);
          for (int v = 0; v < cell_tp::num_vertices_of_face(cell, f); ++v)
            _vtx_to_hfs[cell_vertex(cell, face_vertices[v])] = hf_handle_t(cell_handle, f);
        }
    }
  }

  template<class R, class CT, class OP, class TP>
  std::tuple<R, R, R, R, R, R> parallel_mesh_3d<R, CT, OP, TP>::local_bounding_box() const
  {
    R xmin = std::numeric_limits<R>::max();
    R ymin = std::numeric_limits<R>::max();
    R zmin = std::numeric_limits<R>::max();
    R xmax = std::numeric_limits<R>::lowest();
    R ymax = std::numeric_limits<R>::lowest();
    R zmax = std::numeric_limits<R>::lowest();

    for (size_type c = 0; c < _local_cells.size(); ++c)
    {
      const cell_type& cell = _local_cells[c];
      for (int v = 0; v < cell_tp::num_vertices(cell); ++v)
      {
        const point_type& point = _vertices[cell_vertex(cell, v)];
        if (point.x < xmin) xmin = point.x;
        if (point.x > xmax) xmax = point.x;
        if (point.y < ymin) ymin = point.y;
        if (point.y > ymax) ymax = point.y;
        if (point.z < zmin) zmin = point.z;
        if (point.z > zmax) zmax = point.z;
      }
    }

    return std::make_tuple(xmin, ymin, zmin, xmax, ymax, zmax);
  }

  template<class R, class CT, class OP, class TP>
  std::tuple<R, R, R> parallel_mesh_3d<R, CT, OP, TP>::cell_centroid(cell_handle_t cell_handle) const
  {
    R x{0}, y{0}, z{0};

    const cell_type& cell = get_cell(cell_handle);
    int num_vertices = cell_tp::num_vertices(cell);
    assert(num_vertices > 0);

    for (int v = 0; v < num_vertices; ++v)
    {
      const point_type& point = _vertices[cell_vertex(cell, v)];
      x += point.x;
      y += point.y;
      z += point.z;
    }

    return std::make_tuple(x / num_vertices, y / num_vertices, z / num_vertices);
  }

  template<class R, class CT, class OP, class TP> template<class InputIt>
  void parallel_mesh_3d<R, CT, OP, TP>::repartition(MPI_Comm comm, InputIt first, InputIt last)
  {
    int num_procs, rank;
    MPI_Comm_size(comm, &num_procs);
    MPI_Comm_rank(comm, &rank);

    // prepare scheme and cell indices

//    send_recv_cells(comm, , , , );

    // NOTE: a similar function, send_recv_cell_properties(), could be added and
    // NOTE: called here to handle the re-distribution of cell-based properties

    // replace data structures

    _ghost_cells.clear();
    _vtx_to_hfs.clear();

//    remove_duplicate_vertices(, );

    construct_topology();
  }

  template<class R, class CT, class OP, class TP>
  void parallel_mesh_3d<R, CT, OP, TP>::construct_ghost_layer()
  {

    // NOTE: a similar function, send_recv_cell_properties(), could be added and
    // NOTE: called here to handle the distribution of cell-based properties for
    // NOTE: the ghost layer

    construct_topology();
  }

  // it0: beginning of send scheme; it1: beginning of indices of cells to send;
  // vtx_it: output of received vertices; cell_it: output of received cell connectivities
  template<class R, class CT, class OP, class TP>
  template<typename RandIt, typename OutIt1, typename OutIt2>
  void parallel_mesh_3d<R, CT, OP, TP>::send_recv_cells(MPI_Comm comm, RandIt it0, RandIt it1,
                                                        OutIt1 vtx_it, OutIt2 cell_it)
  {
    int num_procs;
    MPI_Comm_size(comm, &num_procs);

    // 1. send-receive schemes for vertices
    std::vector<integer_type> vtx_recv_scheme(num_procs);
    std::vector<integer_type> vtx_send_scheme;
    std::vector<integer_type> vtx_to_send;

    RandIt cid_begin = it1;
    for (RandIt it = it0; it < it0 + num_procs; ++it)
    {
      integer_type cid_size = *it;

      std::set<integer_type> vtx_indices;
      for (RandIt c = cid_begin; c != cid_begin + cid_size; ++c)
      {
        const cell_type& cell = _local_cells[*c];
        for (int v = 0; v < cell_tp::num_vertices(cell); ++v)
          vtx_indices.insert(cell_vertex(cell, v));
      }

      vtx_send_scheme.push_back(vtx_indices.size());

      size_type cur_size = vtx_to_send.size();
      for (auto sit = vtx_indices.cbegin(); sit != vtx_indices.cend(); ++sit)
        vtx_to_send.push_back(*sit);
      std::sort(vtx_to_send.begin() + cur_size, vtx_to_send.end()); // may do nothing

      cid_begin += cid_size;
    }

    MPI_Alltoall(vtx_send_scheme.data(), 1, mpi_datatype<integer_type>, vtx_recv_scheme.data(),
                 1, mpi_datatype<integer_type>, comm);

    // 2. send-receive vertices, in the form of points (sets of unique vertices are
    //    sent from each rank but on the receiving rank, there could be duplicates)
    integer_type vtx_recv_size = std::reduce(vtx_recv_scheme.cbegin(), vtx_recv_scheme.cend(), 0);
    std::vector<R> vtx_recv_buffer(3 * vtx_recv_size);

    //    receiving
    int num_recvs = 0;
    std::vector<MPI_Request> recv_requests(num_procs);
    integer_type vtx_offset = 0;
    for (int p = 0; p < num_procs; ++p)
    {
      integer_type count = vtx_recv_scheme[p];
      if (count > 0)
      {
        MPI_Irecv(vtx_recv_buffer.data() + vtx_offset, 3 * count, mpi_datatype<R>, p, 0, comm, &recv_requests[num_recvs++]);
        vtx_offset += (3 * count);
      }
    }

    //    sending
    int num_sends = 0;
    std::vector<MPI_Request> send_requests(num_procs);
    vtx_offset = 0;
    for (int p = 0; p < num_procs; ++p)
    {
      integer_type count = vtx_send_scheme[p];
      if (count > 0)
      {
        std::vector<R> vtx_send_buffer;
        for (integer_type i = 0; i < count; ++i)
        {
          const point_type& point = _vertices[vtx_to_send[i + vtx_offset]];
          vtx_send_buffer.push_back(point.x);
          vtx_send_buffer.push_back(point.y);
          vtx_send_buffer.push_back(point.z);
        }

        MPI_Isend(vtx_send_buffer.data(), 3 * count, mpi_datatype<R>, p, 0, comm, &send_requests[num_sends++]);
        vtx_offset += count;
      }
    }

    std::vector<MPI_Status> sends(num_sends);
    MPI_Waitall(num_sends, send_requests.data(), sends.data());
    std::vector<MPI_Status> recvs(num_recvs);
    MPI_Waitall(num_recvs, recv_requests.data(), recvs.data());

    // 3. send-receive schemes for cell connectivities
    std::vector<integer_type> conn_recv_sizes(num_recvs);
    num_recvs = 0;
    for (int p = 0; p < num_procs; ++p)
      if (vtx_recv_scheme[p] > 0)
      {
        MPI_Irecv(conn_recv_sizes.data() + num_recvs, 1, mpi_datatype<integer_type>, p, 0, comm, &recv_requests[num_recvs]);
        ++num_recvs;
      }

    num_sends = 0;
    cid_begin = it1;
    for (int p = 0; p < num_procs; ++p)
      if (vtx_send_scheme[p] > 0)
      {
        integer_type cid_size = *(it0 + p);
        assert(cid_size > 0);

        integer_type conn_size = 0;
        for (RandIt c = cid_begin; c != cid_begin + cid_size; ++c)
        {
          const cell_type& cell = _local_cells[*c];
          conn_size += cell_tp::num_vertices(cell) + 1; // an extra integer to indicate the cell shape
        }

        MPI_Isend(&conn_size, 1, mpi_datatype<integer_type>, p, 0, comm, &send_requests[num_sends++]);

        cid_begin += cid_size;
      }

    MPI_Waitall(num_sends, send_requests.data(), sends.data());
    MPI_Waitall(num_recvs, recv_requests.data(), recvs.data());

    // 4. send-receive cell connectivities
    integer_type conn_recv_size = std::reduce(conn_recv_sizes.cbegin(), conn_recv_sizes.cend(), 0);
    std::vector<integer_type> conn_recv_buffer(conn_recv_size);

    //    receiving
    num_recvs = 0;
    integer_type conn_offset = 0;
    for (int p = 0; p < num_procs; ++p)
      if (vtx_recv_scheme[p] > 0)
      {
        integer_type count = conn_recv_sizes[num_recvs];
        MPI_Irecv(conn_recv_buffer.data() + conn_offset, count, mpi_datatype<integer_type>, p, 0, comm, &recv_requests[num_recvs]);
        ++num_recvs;
        conn_offset += count;
      }

    //    sending
    num_sends = 0;
    size_type sct_begin = 0;
    cid_begin = it1;
    for (int p = 0; p < num_procs; ++p)
    {
      integer_type sct_size = vtx_send_scheme[p];
      assert(sct_size >= 0);

      if (sct_size > 0)
      {
        integer_type cid_size = *(it0 + p);
        assert(cid_size > 0);

        std::vector<integer_type> conn_send_buffer;
        for (RandIt c = cid_begin; c != cid_begin + cid_size; ++c)
        {
          const cell_type& cell = _local_cells[*c];
          conn_send_buffer.push_back(cell_shape(cell)); // shape_3d enum sent as integer
          for (int v = 0; v < cell_tp::num_vertices(cell); ++v)
          {
            // convert the original vertex index to new index to
            // the vector of vertices sent in 2. NOTE that indices
            // in the section [sct_begin, sct_begin + sct_size) of
            // vtx_to_send are already ordered due to the call to
            // std::sort() in 1 so a binary search will work.
            integer_type vid = cell_vertex(cell, v);
            auto it_srch = std::lower_bound(vtx_to_send.begin() + sct_begin,
                                            vtx_to_send.begin() + sct_begin + sct_size, vid);
            assert(it_srch != vtx_to_send.begin() + sct_begin + sct_size);
            conn_send_buffer.push_back(it_srch - vtx_to_send.begin() - sct_begin);
          }
        }

        MPI_Isend(conn_send_buffer.data(), conn_send_buffer.size(), mpi_datatype<integer_type>,
                  p, 0, comm, &send_requests[num_sends++]);

        sct_begin += sct_size;
        cid_begin += cid_size;
      }
    }

    MPI_Waitall(num_sends, send_requests.data(), sends.data());
    MPI_Waitall(num_recvs, recv_requests.data(), recvs.data());

    // 5. populate output data structures
    for (integer_type i = 0; i < 3 * vtx_recv_size; i += 3)
      *vtx_it = point_3d(vtx_recv_buffer[i], vtx_recv_buffer[i + 1],  vtx_recv_buffer[i + 2]); 

    num_recvs = 0;
    vtx_offset = 0;
    conn_offset = 0;
    for (int p = 0; p < num_procs; ++p)
    {
      integer_type vtx_count = vtx_recv_scheme[p];
      if (vtx_count > 0)
      {
        integer_type conn_count = conn_recv_sizes[num_recvs++];
        for (integer_type i = conn_offset; i < (conn_offset + conn_count); )
        {
          cell_type cell = cell_op::create_cell(static_cast<shape_3d>(conn_recv_buffer[i++]));
          for (int v = 0; v < cell_tp::num_vertices(cell); ++v)
            cell_vertex(cell, v, vtx_offset + conn_recv_buffer[i++]);
          *cell_it = cell;
        }

        vtx_offset += vtx_count;
        conn_offset += conn_count;
      }
    }
  }

  template<class R, class CT, class OP, class TP>
  void parallel_mesh_3d<R, CT, OP, TP>::remove_duplicate_vertices(std::vector<point_type>& vertices,
                                                                  std::vector<cell_type>& cells)
  {
    using vertex_comparer = vertex_3d_comparer<R, integer_type>; // only compare x, y, and z, not id
    using vertex_type = vertex_3d<R, integer_type>;
    kd_tree<3, vertex_type, vertex_comparer> vtx_tree;

    std::vector<std::pair<integer_type, integer_type>> vtx_to_remove;
    for (size_type i = 0; i < vertices.size(); ++i)
    {
      const point_type& point = vertices[i];
      vertex_type vertex{point.x, point.y, point.z, i};
      auto it = vtx_tree.find(vertex);
      if (it == vtx_tree.end()) vtx_tree.insert(vertex);
      else vtx_to_remove.push_back(std::make_pair(i, it->id));
    }

    integer_type num_removed = 0;
    for (auto vit = vtx_to_remove.cbegin(); vit != vtx_to_remove.cend(); ++vit)
    {
      integer_type removed = vit->first;
      integer_type replacement = vit->second;

      vertices.erase(vertices.begin() + removed - num_removed);
      ++num_removed;

      for (auto cit = cells.begin(); cit != cells.end(); ++cit)
      {
        cell_type& cell = *cit;
        for (int v = 0; v < cell_tp::num_vertices(cell); ++v)
        {
          integer_type vid = cell_vertex(cell, v);
          if (vid == removed) cell_vertex(cell, v, replacement);
          else if (vid > removed) cell_vertex(cell, v, vid - 1);
        }
      }
    }
  }

  template<class R, class CT, class OP, class TP>
  void parallel_mesh_3d<R, CT, OP, TP>::print_connectivity(std::ostream& out) const
  {
    out << "local cell-to-vertex table (size = " << _local_cells.size() << ")" << std::endl;
    for (size_type c = 0; c < _local_cells.size(); ++c)
    {
      const cell_type& cell = _local_cells[c];
      out << c << ": ";
      for (int v = 0; v < cell_tp::num_vertices(cell); ++v)
        out << "   " << cell_vertex(cell, v);
      out << std::endl;
    }

    out << "ghost cell-to-vertex table (size = " << _ghost_cells.size() << ")" << std::endl;
    for (size_type c = 0; c < _ghost_cells.size(); ++c)
    {
      const cell_type& cell = _ghost_cells[c];
      out << c << ": ";
      for (int v = 0; v < cell_tp::num_vertices(cell); ++v)
        out << "   " << cell_vertex(cell, v);
      out << std::endl;
    }
  }

  template<class R, class CT, class OP, class TP>
  void parallel_mesh_3d<R, CT, OP, TP>::print_twin_hfs(std::ostream& out) const
  {
    out << "local twin-half-facet table (size = " << _local_cells.size() << ")" << std::endl;
    for (size_type c = 0; c < _local_cells.size(); ++c)
    {
      const cell_type& cell = _local_cells[c];
      out << c << ":";
      for (int f = 0; f < cell_tp::num_faces(cell); ++f)
      {
        hf_handle_t hf = twin_hf(cell, f);
        if (hf.is_valid())
        {
          if (hf.cell_handle().is_ghost()) out << "  <g";
          else out << "  <";
          out << hf.cell_handle().local_index() << ", " << hf.face_index() << ">";
        }
        else out << "  <-->";
      }
      out << std::endl;
    }

    out << "ghost twin-half-facet table (size = " << _ghost_cells.size() << ")" << std::endl;
    for (size_type c = 0; c < _ghost_cells.size(); ++c)
    {
      const cell_type& cell = _ghost_cells[c];
      out << c << ":";
      for (int f = 0; f < cell_tp::num_faces(cell); ++f)
      {
        hf_handle_t hf = twin_hf(cell, f);
        if (hf.is_valid())
        {
          if (hf.cell_handle().is_ghost()) out << "  <g";
          else out << "  <";
          out << hf.cell_handle().local_index() << ", " << hf.face_index() << ">";
        }
        else out << "  <-->";
      }
      out << std::endl;
    }
  }

  template<class R, class CT, class OP, class TP>
  void parallel_mesh_3d<R, CT, OP, TP>::export_to_msh_format(std::ostream& out, int config) const
  {
    assert(config >= 0 && config <= 2);

    out << "$MeshFormat" << std::endl;
    out << "4.1 0 8" << std::endl;
    out << "$EndMeshFormat" << std::endl;
	
    // Nodes
    out << "$Nodes" << std::endl;
    out << "1 " << num_vertices() << " 1 " << num_vertices() << std::endl;
    out << "3 1 0 " << num_vertices() << std::endl;
    for (size_type v = 1; v <= _vertices.size(); ++v) // 1-based index in Gmsh!
      out << v << std::endl;
    for (size_type v = 1; v <= _vertices.size(); ++v)
    {
      const point_type& point = _vertices[v - 1]; // 1-based index in Gmsh!
      out << point.x << " " << point.y << " " << point.z << std::endl;
    }
    out << "$EndNodes" << std::endl;

    // Elements - TODO: group cells based on shape and write each group
    //            TODO: in a different section
    out << "$Elements" << std::endl;
    int num_cells = config == 0 ? num_local_cells() : (config == 1 ? num_ghost_cells() :
                    num_local_cells() + num_ghost_cells());
    out << "1 " << num_cells << " 1 " << num_cells << std::endl;
    out << "3 1 5 " << num_cells << std::endl; // NOTE: HARD-CODED SHAPE OF HEX (ELEMENT TYPE = 5)!
        
    int cell_index = 1; // 1-based index in Gmsh!
    if (config == 0 || config == 2)
      for (auto it = _local_cells.cbegin(); it != _local_cells.cend(); ++it)
      {
        const cell_type& cell = *it;
        out << cell_index++;
        for (int v = 0; v < cell_tp::num_vertices(cell); ++v)
                out << " " << cell_vertex(cell, v) + 1; // vertex index is also 1-based!
        out << std::endl;
      }

    if (config == 1 || config == 2)
      for (auto it = _ghost_cells.cbegin(); it != _ghost_cells.cend(); ++it)
      {
        const cell_type& cell = *it;
        out << cell_index++;
        for (int v = 0; v < cell_tp::num_vertices(cell); ++v)
                out << " " << cell_vertex(cell, v) + 1; // vertex index is also 1-based!
        out << std::endl;
      }

    out << "$EndElements" << std::endl;
  }

} // namespace pmh

#endif
