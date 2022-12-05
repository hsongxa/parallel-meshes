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
#include "kd_tree.h"
#include "cell_types.h"
#include "cell_type_traits.h"

#include <cstddef>
#include <vector>
#include <map>
#include <cassert>
#include <algorithm>
#include <iostream>

namespace pmh {

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
  // construction of ghost layers is performed after the
  // partitioning/repartitioning.
  //
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

    explicit parallel_mesh_3d(int part_id) : _part_id(part_id) {}

    // ------ local construction operations ------

    // population of vertices used by local cells (no ghosts)
    template<class InputIt>
    void fill_local_vertices(InputIt first, InputIt last);

    // population of local cell connectivities
    template<class InputIt>
    void fill_local_cells(InputIt first, InputIt last);

    // ------ MPI operations ------

    // partitioning/repartitioning
    template<class InputIt>
    void repartition(InputIt first, InputIt last);

    // one layer of ghost cells that are face neighbors of local cells
    void construct_ghost_layer();

    // ------ queries of vertices, faces, and cells ------

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

    // ------ helpers for mesh I/O ------

    // The following two functions are the ONLY places where we must be explicit
    // about the cell types, which are needed when doing mesh I/O
    static cell_type& emplace_default_cell(std::vector<cell_type>& cell_vector, shape_3d shape)
    {
      if (shape == shape_3d::HEXAHEDRON)
        return cell_vector.emplace_back(hex_t());
      else if (shape == shape_3d::WEDGE)
        return cell_vector.emplace_back(wdg_t());
      else if (shape == shape_3d::PYRAMID)
        return cell_vector.emplace_back(prm_t());
      else
      {
        assert(shape == shape_3d::TETRAHEDRON);
        return cell_vector.emplace_back(tet_t());
      }
    }

  private:
    std::vector<point_type>  _vertices;
    std::vector<hf_handle_t> _vtx_to_hfs;

    std::vector<cell_type>   _local_cells;
    std::vector<cell_type>   _ghost_cells;

    int _part_id; // rank

    // only used in MPI operations
    using vertex_comparer = vertex_3d_comparer<R, std::int64_t>;
    kd_tree<3, vertex_3d<R, std::int64_t>, vertex_comparer> _vtx_tree;

  public: // mainly for debugging
    void print_connectivity(std::ostream& out) const;
    void print_twin_hfs(std::ostream& out) const;
    void export_to_msh_format(std::ostream& out, int config) const;
  };

  template<class R, class CT, class OP, class TP> template<class InputIt>
  void parallel_mesh_3d<R, CT, OP, TP>::fill_local_vertices(InputIt first, InputIt last)
  {
    _vtx_to_hfs.clear();
    _vertices.clear();

    while (first != last)
      _vertices.emplace_back(point_3d{*first++, *first++, *first++});
  }

  template<class R, class CT, class OP, class TP> template<class InputIt>
  void parallel_mesh_3d<R, CT, OP, TP>::fill_local_cells(InputIt first, InputIt last)
  {
    _ghost_cells.clear();
    _local_cells.clear();

    while (first != last) // format of each cell is: <shape> <v1, v2, ...>
    {
      shape_3d shape = static_cast<shape_3d>(*first++);
      cell_type& cell = emplace_default_cell(_local_cells, shape);
      for (int i = 0; i < cell_tp::num_vertices(cell); ++i)
        cell_vertex(cell, i, *first++);
    }
  }

  template<class R, class CT, class OP, class TP>
  void parallel_mesh_3d<R, CT, OP, TP>::fill_twin_hfs()
  {
    int face_vertices[4];  // maximum 4 vertices in each face for all cell types
    std::map<vid_tuple<integer_type>, hf_handle_t> dict;

    for (std::size_t c = 0; c < _local_cells.size(); ++c)
    {
      cell_type& cell = _local_cells[c];
      cell_handle_t cell_handle = cell_handle_t(c, false);
      for (int f = 0; f < cell_tp::num_faces(cell); ++f)
      {
        cell_tp::vertices_of_face(cell, f, face_vertices);
        int num_vertices = cell_tp::num_vertices_of_face(cell, f);

        vid_tuple<integer_type> vtuple = num_vertices < 4 ?
          std::make_tuple(-1, cell_vertex(cell, face_vertices[0]),
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
          dict.insert(std::make_pair(vtuple, hf_handle_t(cell_handle, f)));
      }
    }

    for (std::size_t c = 0; c < _ghost_cells.size(); ++c)
    {
      cell_type& cell = _ghost_cells[c];
      cell_handle_t cell_handle = cell_handle_t(c, true);
      for (int f = 0; f < cell_tp::num_faces(cell); ++f)
      {
        cell_tp::vertices_of_face(cell, f, face_vertices);
        int num_vertices = cell_tp::num_vertices_of_face(cell, f);

        vid_tuple<integer_type> vtuple = num_vertices < 4 ?
          std::make_tuple(-1, cell_vertex(cell, face_vertices[0]),
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
          dict.insert(std::make_pair(vtuple, hf_handle_t(cell_handle, f)));
      }
    }
  }

  // ASSUMPTION: this is always called after fill_twin_hfs()
  template<class R, class CT, class OP, class TP>
  void parallel_mesh_3d<R, CT, OP, TP>::fill_vtx_to_hfs()
  {
    int incident_faces[4]; // maximum 4 faces meet at a vertex for all cell types
    int face_vertices[4];  // maximum 4 vertices in each face for all cell types

    // allocate and initialize with invalid hfs
    _vtx_to_hfs.resize(_vertices.size());
    std::fill(_vtx_to_hfs.begin(), _vtx_to_hfs.end(), hf_handle_t());

    // first, go through ghost cells
    for (std::size_t c = 0; c < _ghost_cells.size(); ++c)
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
    for (std::size_t c = 0; c < _local_cells.size(); ++c)
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

  template<class R, class CT, class OP, class TP> template<class InputIt>
  void parallel_mesh_3d<R, CT, OP, TP>::repartition(InputIt first, InputIt last)
  {
  }

  template<class R, class CT, class OP, class TP>
  void parallel_mesh_3d<R, CT, OP, TP>::construct_ghost_layer()
  {
  }

  template<class R, class CT, class OP, class TP>
  void parallel_mesh_3d<R, CT, OP, TP>::print_connectivity(std::ostream& out) const
  {
  }

  template<class R, class CT, class OP, class TP>
  void parallel_mesh_3d<R, CT, OP, TP>::print_twin_hfs(std::ostream& out) const
  {
  }

  template<class R, class CT, class OP, class TP>
  void parallel_mesh_3d<R, CT, OP, TP>::export_to_msh_format(std::ostream& out, int config) const
  {
  }

/*
template<class CT0, class CT1, class CT2, class CT3>
void parallel_mesh<CT0, CT1, CT2, CT3>::print_connectivity_table(std::ostream& out) const
{
	out << "local cell-to-vertex table (size = " << _local_cells.size() << ", used size = " << _local_cells.used_size() << ")" << std::endl;
	for (auto it = _local_cells.cbegin(); it != _local_cells.cend(); ++it)
	{
		const cell_type& cell = *it;
		out << it.at_index() << ": * " << cell.index() << " *";
		for (int v = 0; v < reference_shape_3d::num_vertices(cell_shape(cell)); ++v)
			out << "   " << cell_vertex(cell, v);
		out << std::endl;
	}
	out << "ghost cell-to-vertex table (size = " << _ghost_cells.size() << ", used size = " << _ghost_cells.used_size() << ")" << std::endl;
	for (auto it = _ghost_cells.cbegin(); it != _ghost_cells.cend(); ++it)
	{
		const cell_type& cell = *it;
		out << it.at_index() << ": * " << cell.index() << " *";
		for (int v = 0; v < reference_shape_3d::num_vertices(cell_shape(cell)); ++v)
			out << "   " << cell_vertex(cell, v);
		out << std::endl;
	}
}

template<class CT0, class CT1, class CT2, class CT3>
void parallel_mesh<CT0, CT1, CT2, CT3>::print_twin_hfs_table(std::ostream& out) const
{
	out << "local twin-half-facet table (size = " << _local_cells.size() << ", used size = " << _local_cells.used_size() << ")" << std::endl;
	for (auto it = _local_cells.cbegin(); it != _local_cells.cend(); ++it)
	{
		const cell_type& cell = *it;
		out << it.at_index() << ":";
		for (int f = 0; f < reference_shape_3d::num_faces(cell_shape(cell)); ++f)
		{
			hf_handle_t hf = twin_hf(cell, f);
			if (hf.is_valid())
			{
				if (hf.cell_index().is_ghost())
					out << "  <g";
				else
					out << "  <";
				out << hf.cell_index().index() << ", " << hf.face_index() << ">";
			}
			else
				out << "  <-->";
		}
		out << std::endl;
	}
	out << "ghost twin-half-facet table (size = " << _ghost_cells.size() << ", used size = " << _ghost_cells.used_size() << ")" << std::endl;
	for (auto it = _ghost_cells.cbegin(); it != _ghost_cells.cend(); ++it)
	{
		const cell_type& cell = *it;
		out << it.at_index() << ":";
		for (int f = 0; f < reference_shape_3d::num_faces(cell_shape(cell)); ++f)
		{
			hf_handle_t hf = twin_hf(cell, f);
			if (hf.is_valid())
			{
				if (hf.cell_index().is_ghost())
					out << "  <g";
				else
					out << "  <";
				out << hf.cell_index().index() << ", " << hf.face_index() << ">";
			}
			else
				out << "  <-->";
		}
		out << std::endl;
	}
}

template<class CT0, class CT1, class CT2, class CT3>
void parallel_mesh<CT0, CT1, CT2, CT3>::export_to_msh_format(std::ostream& out, int config) const
{
	assert(config >= 0 && config <= 2);
	int cell_index = 1; // 1-based index in Gmsh!
	if (config == 0 || config == 2)
	{
		for (auto it = _local_cells.cbegin(); it != _local_cells.cend(); ++it)
		{
			const cell_type& cell = *it;
			out << cell_index++;
			for (int v = 0; v < reference_shape_3d::num_vertices(cell_shape(cell)); ++v)
				out << " " << cell_vertex(cell, v) + 1; // vertex index is also 1-based!
			out << std::endl;
		}
	}
	if (config == 1 || config == 2)
	{
		for (auto it = _ghost_cells.cbegin(); it != _ghost_cells.cend(); ++it)
		{
			const cell_type& cell = *it;
			out << cell_index++;
			for (int v = 0; v < reference_shape_3d::num_vertices(cell_shape(cell)); ++v)
				out << " " << cell_vertex(cell, v) + 1; // vertex index is also 1-based!
			out << std::endl;
		}
	}
}
*/
} // namespace pmh

#endif
