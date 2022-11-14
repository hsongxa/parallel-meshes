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

#include "vertex_type.h"
#include "cell_types.h"
#include "cell_operators.h"

#include <cstddef>
#include <vector>
#include <cassert>
#include <iostream>

namespace pmh {

  // R - floating point type for coordinates, CT - cell type,
  // OP - cell operator type
  //
  // Passing a tet_t type to CT will make a tetrahedral mesh,
  // a hex_t for hexahedral mesh, ..., etc. For mixed
  // tet-wdg-hex mesh, CT is std::variant<hex_t, wdg_t, tet_t>.
  //
  // Similarly, std::variant<hex_t, wdg_t, prm_t, tet_t> will
  // give a mesh with mixed cell shapes of all four.
  template< class R, class CT,
            class OP = typename cell_operator_traits<CT>::operator_type >
  class parallel_mesh_3d
  {
  public:
    using vertex_type = vertex_3d<R, std::int64_t>;
    using cell_type = CT;
    using cell_op = OP;
    using size_type = std::size_t;

    explicit parallel_mesh_3d(int part_id) : _part_id(part_id) {}

    // ------ local construction operations ------

    // population of vertices used by local cells (no ghosts)
    template<class InputIt>
    void fill_vertices(InputIt first, InputIt last)
    {
      integer_type index = 0;
      while (first != last)
        _vertices.emplace_back(vertex_type{*first++, *first++, *first++, index++});
    }

    // population of local cell connectivities
    template<class InputIt>
    void fill_connectivity(InputIt first, InputIt last);

    // population of twin half-facets among local cells
    // NOTE: twin half-facets between local and ghost cells
    // NOTE: are not populated by this funciton -- they are
    // NOTE: the results of construct_ghost_layer();
    // NOTE: (twin half-facets among ghost cells themselves
    // NOTE: are never populated in this implementation)
    void fill_local_twin_hfs();

    // ------ MPI operations ------

    // partitioning/repartitioning
    template<class InputIt>
    void repartition(InputIt first, InputIt last);

    // one layer of ghost cells that are face neighbors of local cells
    void construct_ghost_layer();

    // ------ queries ------

    size_type num_vertices() const { return _vertices.size(); }

    size_type num_local_cells() const { return _local_cells.size(); }

    size_type num_ghost_cells() const { return _ghost_cells.size(); }

    const vertex_type& get_vertex(integer_type vid) const { return _vertices[vid]; }

    vertex_type& get_vertex(integer_type vid) { return _vertices[vid]; }

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

    // invalid cell handle is returned for boundary faces
    cell_handle_t get_face_neighbor(cell_handle_t cell_handle, int f) const
    {
      hf_handle_t hf = cell_op::twin_hf(get_cell(cell_handle), f);
      if (!hf.is_valid()) return cell_handle_t();
      return hf.cell_handle();
    }

    // ------ static queries on a given cell ------

    static shape_3d cell_shape(const cell_type& cell) { return cell_op::shape(cell); }

    static integer_type cell_vertex(const cell_type& cell, int v)
      { return cell_op::cell_vertex(cell, v); }

    static hf_handle_t twin_hf(const cell_type& cell, int f)
      { return cell_op::twin_hf(cell, f); }

  private:
    static void cell_vertex(cell_type& cell, int v, integer_type vertex_id)
      { cell_op::cell_vertex(cell, v, vertex_id); }

    static void twin_hf(cell_type& cell, int f, hf_handle_t hf)
      { cell_op::twin_hf(cell, f, hf); }

/*
	// for triangular face, the first item of the returned tuple is -1
	static util::index_tuple vertices_of_face(const cell_type& cell, shape_3d shape, int face, int face_vertices[])
	{
		reference_shape_3d::vertices_of_face(shape, face, face_vertices);
		if (reference_shape_3d::num_vertices_of_face(shape, face) == 3)
			return std::make_tuple(-1, cell_vertex(cell, face_vertices[0]),	cell_vertex(cell, face_vertices[1]),
								   cell_vertex(cell, face_vertices[2]));
		return std::make_tuple(cell_vertex(cell, face_vertices[0]), cell_vertex(cell, face_vertices[1]),
							   cell_vertex(cell, face_vertices[2]), cell_vertex(cell, face_vertices[3]));
	}
*/
  private:
    std::vector<cell_type>   _local_cells;
    std::vector<cell_type>   _ghost_cells;
    std::vector<vertex_type> _vertices;

    int _part_id; // rank

  public: // mainly for debugging
    void print_connectivity_table(std::ostream& out) const;
    void print_twin_hfs_table(std::ostream& out) const;
    void export_to_msh_format(std::ostream& out, int config) const;
  };


  template<class R, class CT, class OP> template<class InputIt>
  void parallel_mesh_3d<R, CT, OP>::fill_connectivity(InputIt first, InputIt last)
  {
  }

  template<class R, class CT, class OP>
  void parallel_mesh_3d<R, CT, OP>::fill_local_twin_hfs()
  {
  }

  template<class R, class CT, class OP> template<class InputIt>
  void parallel_mesh_3d<R, CT, OP>::repartition(InputIt first, InputIt last)
  {
  }

  template<class R, class CT, class OP>
  void parallel_mesh_3d<R, CT, OP>::construct_ghost_layer()
  {
  }

  template<class R, class CT, class OP>
  void parallel_mesh_3d<R, CT, OP>::print_connectivity_table(std::ostream& out) const
  {
  }

  template<class R, class CT, class OP>
  void parallel_mesh_3d<R, CT, OP>::print_twin_hfs_table(std::ostream& out) const
  {
  }

  template<class R, class CT, class OP>
  void parallel_mesh_3d<R, CT, OP>::export_to_msh_format(std::ostream& out, int config) const
  {
  }
/*
template<class CT0, class CT1, class CT2, class CT3> template<class InputIt>
void parallel_mesh<CT0, CT1, CT2, CT3>::fill_connectivity(InputIt first, InputIt last, bool are_ghost)
{
	assert(_level == 0); // only base level needs to do this
	discontinuous_vector<cell_type>& cell_vec = are_ghost ? _ghost_cells : _local_cells;
	cell_vec.clear();

	while (first != last)
	{
		shape_3d shape = static_cast<shape_3d>(*first++);
		cell_type& cell = emplace_default_cell(cell_vec, shape);
		if (are_ghost) tree_id(cell, *first++); // directly fill ghost cells' global id
		for (int i = 0; i < reference_shape_3d::num_vertices(shape); ++i)
			cell_vertex(cell, i, *first++);
	}
}

template<class CT0, class CT1, class CT2, class CT3>
void parallel_mesh<CT0, CT1, CT2, CT3>::fill_local_twin_hfs()
{
	assert(_level == 0); // only base level needs to do this

	// NOTE that this is based on the assumption that the base grid
	// is conformal, i.e., there are two and only two cells sharing
	// a common face except boundary faces. If this is not true,
	// then a multimap should be used and the following code modified
	std::map<util::index_tuple, hf_handle_t> dict;

	int local_vertices[reference_shape_3d::MAX_NUM_VERTICES_PER_FACE];
	integer_type vtx_indices[reference_shape_3d::MAX_NUM_VERTICES_PER_FACE];
	util::index_tuple ind_tuple;
	for (auto it = _local_cells.begin(); it != _local_cells.end(); ++it)
	{
		cell_type& cell = *it;
		shape_3d shape = cell_shape(cell);
		for (int f = 0; f < reference_shape_3d::num_faces(shape); ++f)
		{
			reference_shape_3d::vertices_of_face(shape, f, local_vertices);
			int num = reference_shape_3d::num_vertices_of_face(shape, f);
			for (int i = 0; i < num; ++i)
				vtx_indices[i] = cell_vertex(cell, local_vertices[i]);

			if (num < reference_shape_3d::MAX_NUM_VERTICES_PER_FACE)
			{
				assert(num == 3);
				ind_tuple = std::make_tuple(-1, vtx_indices[0], vtx_indices[1], vtx_indices[2]);
			}
			else
				ind_tuple = std::make_tuple(vtx_indices[0], vtx_indices[1], vtx_indices[2], vtx_indices[3]);
			util::order_tuple_items(ind_tuple);

			// if already exists a hf with the same key, then they are twin hfs
			auto it_twin = dict.find(ind_tuple);
			if (it_twin != dict.end())
			{
				hf_handle_t hf = it_twin->second;
				twin_hf(cell, f, hf);
				twin_hf(get_cell(hf.cell_index()),  hf.face_index(), hf_handle_t(_level, it.at_index(), false, f));
				dict.erase(it_twin);
			}
			// otherwise
			else
				dict.insert(std::make_pair(ind_tuple, hf_handle_t(_level, it.at_index(), false, f)));
		}
	}

	// the remaining "single" facets need to be sent to other ranks
	// to find their "pairs", hence the ghost cells...

	// after the above, the still remaining "single" facets are boundary facets...
}

template<class CT0, class CT1, class CT2, class CT3>
bool parallel_mesh<CT0, CT1, CT2, CT3>::is_twin_edge_refined(const cell_type& cell, int local_edge, const cell_type** twin_cell, int* twin_edge) const
{
	// get the global indices of the vertices defining the edge
	shape_3d shape = cell_shape(cell);
	auto local_vertices = reference_shape_3d::vertices_of_edge(shape, local_edge);
	const auto vertex0 = cell_vertex(cell, local_vertices.first);
	const auto vertex1 = cell_vertex(cell, local_vertices.second);
	
	// NOTE + TODO: NEED TO SEARCH IN BOTH DIRECTIONS BECAUSE BOUNDARY MAY HAPPEN IN THE MIDDLE OF THE SEARCH!!!

	// the direction along which to visit neighboring cells is chosen to avoid boundary facet
	auto local_faces = reference_shape_3d::incident_faces(shape, local_edge);
	hf_handle_t nxt_hf_handle = twin_hf(cell, local_faces.first);
	if (!nxt_hf_handle.is_valid()) // boundary facet
		nxt_hf_handle = twin_hf(cell, local_faces.second);

	// visit neighboring cells incident to the edge
	int face_edges[reference_shape_3d::MAX_NUM_VERTICES_PER_FACE]; // number of face edges is the same as number of face vertices
	while (nxt_hf_handle.is_valid())
	{
		const cell_type& nxt_cell = get_cell(nxt_hf_handle.cell_index());
		if (&cell == &nxt_cell) break;

		// search for the local edge in the nxt_cell's nxt_face, whose end vertices are [vertex0, vertex1]
		auto nxt_shape = cell_shape(nxt_cell);
		int nxt_face = nxt_hf_handle.face_index();
		int nxt_edge = -1;
		reference_shape_3d::edges_of_face(nxt_shape, nxt_face, face_edges);
		for (int i = 0; i < reference_shape_3d::num_edges_of_face(nxt_shape, nxt_face); ++i)
		{
			local_vertices = reference_shape_3d::vertices_of_edge(nxt_shape, face_edges[i]);
			if ((cell_vertex(nxt_cell, local_vertices.first) == vertex0 && cell_vertex(nxt_cell, local_vertices.second) == vertex1) ||
				(cell_vertex(nxt_cell, local_vertices.first) == vertex1 && cell_vertex(nxt_cell, local_vertices.second) == vertex0))
			{
				nxt_edge = face_edges[i];
				break;
			}
		}
		assert(nxt_edge >= 0);

		if (is_cell_refined(nxt_cell))
		{
			*twin_cell = &nxt_cell;
			*twin_edge = nxt_edge;
			return true;
		}

		local_faces = reference_shape_3d::incident_faces(nxt_shape, nxt_edge);
		assert(local_faces.first == nxt_face || local_faces.second == nxt_face);
		nxt_hf_handle = twin_hf(nxt_cell, local_faces.first == nxt_face ? local_faces.second : local_faces.first);
	}

	return false;
}

template<class CT0, class CT1, class CT2, class CT3>
bool parallel_mesh<CT0, CT1, CT2, CT3>::is_twin_hf_refined(const cell_type& cell, int local_face, const cell_type** twin_cell, int* twin_face) const
{
	hf_handle_t twin_hf_handle = twin_hf(cell, local_face);
	if (twin_hf_handle.is_valid())
	{
		const auto& cell = get_cell(twin_hf_handle.cell_index());
		if (is_cell_refined(cell))
		{
			*twin_cell = &cell;
			*twin_face = twin_hf_handle.face_index();
			return true;
		}
	}
	return false;
}

template<class CT0, class CT1, class CT2, class CT3>
std::pair<cell_handle_t, int> parallel_mesh<CT0, CT1, CT2, CT3>::add_refinement(shape_3d shape, int variation, bool is_ghost, in_tree_id_t parent_in_tree_id, const integer_type vtx_indices[],
																				   int num_faces_to_map, const int faces_to_map[], const cell_type* nb_cells_to_map[], const int nb_faces_to_map[])
{
	assert(_level > 0); // refinement is never added to the base level

	const int num_children = refinement_template_3d::num_children(shape);
	auto& cell_vec = is_ghost ? _ghost_cells : _local_cells;
	const auto first_child_index = cell_vec.require_chunk(num_children);

	// create children cells and fill topology from refinement template
	int local_vertices[reference_shape_3d::MAX_NUM_VERTICES_PER_CELL];
	for (int c = 0; c < num_children; ++c)
	{
		const shape_3d child_shape = refinement_template_3d::child_shape(shape, c);
		cell_type& child = emplace_default_cell(cell_vec, child_shape, c);

		// in_tree_id
		this->in_tree_id(child, parent_in_tree_id.append_refinement(_level, c));

		// cell connectivity
		refinement_template_3d::child_connectivity(shape, c, local_vertices, variation);
		for (int v = 0; v < reference_shape_3d::num_vertices(child_shape); ++v)
			cell_vertex(child, v, vtx_indices[local_vertices[v]]);

		// twin half-facets
		for (int f = 0; f < reference_shape_3d::num_faces(child_shape); ++f)
		{
			const auto pair = refinement_template_3d::child_twin_face(shape, c, f, variation);
			if (pair.first != c) twin_hf(child, f, hf_handle_t(_level, first_child_index + pair.first, is_ghost, pair.second));
		}
	}

	// map children facets with outside cells that were refined
	int children[refinement_template_3d::MAX_NUM_CHILDREN_FACES_PER_FACE], nb_children[refinement_template_3d::MAX_NUM_CHILDREN_FACES_PER_FACE];
	int children_faces[refinement_template_3d::MAX_NUM_CHILDREN_FACES_PER_FACE], nb_children_faces[refinement_template_3d::MAX_NUM_CHILDREN_FACES_PER_FACE];
	int face_vertices[reference_shape_3d::MAX_NUM_VERTICES_PER_FACE];
	int num_children_faces, num_nb_children_faces;
	for (int f = 0; f < num_faces_to_map; ++f)
	{
		const cell_type& nb_cell = *nb_cells_to_map[f];
		const shape_3d nb_shape = cell_shape(nb_cell);
		const int nb_variation = nb_shape == shape_3d::TETRAHEDRON ? std::get<CT3>(nb_cell)._variation : 0;
		assert(nb_variation >= 0 && nb_variation <= 2);

		num_children_faces = refinement_template_3d::children_faces_incident_to_face(shape, faces_to_map[f], children, children_faces, variation);
		num_nb_children_faces = refinement_template_3d::children_faces_incident_to_face(nb_shape, nb_faces_to_map[f], nb_children, nb_children_faces, nb_variation);
		assert(num_children_faces == num_nb_children_faces);

		// match children hfs
		auto anchor = vertices_of_face(cell_vec[first_child_index + children[0]], refinement_template_3d::child_shape(shape, children[0]), children_faces[0], face_vertices);
		util::order_tuple_items(anchor);

		int pos1 = -1;
		for (int i = 0; i < num_nb_children_faces; ++i)
		{
			auto ftuple = vertices_of_face(get_cell(first_child(nb_cell).forward(nb_children[i])),
								refinement_template_3d::child_shape(nb_shape, nb_children[i]), nb_children_faces[i], face_vertices);
			util::order_tuple_items(ftuple);
			if (ftuple == anchor) {	pos1 = i; break; }
		}
		assert(pos1 >= 0);

		if (reference_shape_3d::num_vertices_of_face(shape, faces_to_map[f]) == 3) // triangle hfs
		{
			assert(pos1 < num_children_faces - 1); // the last facet is in the middle - should automatically match
			cell_handle_t nb_child_idx;
			for (int pos0 = 0; pos0 < num_children_faces - 1; ++pos0)
			{
				nb_child_idx = first_child(nb_cell).forward(nb_children[pos1]);
				twin_hf(cell_vec[first_child_index + children[pos0]], children_faces[pos0], hf_handle_t(nb_child_idx, nb_children_faces[pos1]));
				twin_hf(get_cell(nb_child_idx), nb_children_faces[pos1], hf_handle_t(_level, first_child_index + children[pos0], is_ghost, children_faces[pos0]));

				--pos1;
				if (pos1 < 0) pos1 = num_children_faces - 2;
			}
			pos1 = num_children_faces - 1;
			nb_child_idx = first_child(nb_cell).forward(nb_children[pos1]);
			twin_hf(cell_vec[first_child_index + children[pos1]], children_faces[pos1], hf_handle_t(nb_child_idx, nb_children_faces[pos1]));
			twin_hf(get_cell(nb_child_idx), nb_children_faces[pos1], hf_handle_t(_level, first_child_index + children[pos1], is_ghost, children_faces[pos1]));
		}
		else // quad hfs
		{
			assert(reference_shape_3d::num_vertices_of_face(shape, faces_to_map[f]) == 4);
			for (int pos0 = 0; pos0 < num_children_faces; ++pos0)
			{
				cell_handle_t nb_child_idx = first_child(nb_cell).forward(nb_children[pos1]);
				twin_hf(cell_vec[first_child_index + children[pos0]], children_faces[pos0],	hf_handle_t(nb_child_idx, nb_children_faces[pos1]));
				twin_hf(get_cell(nb_child_idx), nb_children_faces[pos1], hf_handle_t(_level, first_child_index + children[pos0], is_ghost, children_faces[pos0]));

				--pos1;
				if (pos1 < 0) pos1 = num_children_faces - 1;
			}
		}
	}

	return std::make_pair(cell_handle_t(first_child_index, _level, is_ghost), num_children);
}

template<class CT0, class CT1, class CT2, class CT3>
void parallel_mesh<CT0, CT1, CT2, CT3>::remove_refinement(cell_handle_t begin_cell_handle, int num_cells,
															   int num_faces_to_unmap, const cell_type* nb_cells_to_unmap[], const int nb_faces_to_unmap[])
{
	assert(begin_cell_handle.is_valid());
	assert(_level > 0 && _level == begin_cell_handle.level()); // refinement was never added to the base level

	// remove cells
	if (begin_cell_handle.is_ghost())
		_ghost_cells.erase_chunk(begin_cell_handle.index(), num_cells);
	else
		_local_cells.erase_chunk(begin_cell_handle.index(), num_cells);

	// unmap facets
	int nb_children[refinement_template_3d::MAX_NUM_CHILDREN_FACES_PER_FACE];
	int nb_children_faces[refinement_template_3d::MAX_NUM_CHILDREN_FACES_PER_FACE];
	for (int f = 0; f < num_faces_to_unmap; ++f)
	{
		const cell_type& nb_cell = *nb_cells_to_unmap[f];
		const shape_3d nb_shape = cell_shape(nb_cell);
		const int nb_variation = nb_shape == shape_3d::TETRAHEDRON ? std::get<CT3>(nb_cell)._variation : 0;
		assert(nb_variation >= 0 && nb_variation <= 2);

		const auto nb_first_child = first_child(nb_cell);
		const int nb_num_children = refinement_template_3d::children_faces_incident_to_face(nb_shape, nb_faces_to_unmap[f],
										nb_children, nb_children_faces, nb_variation);
		assert(nb_first_child.is_valid() && nb_num_children > 0);
		for (int i = 0; i < nb_num_children; ++i)
			twin_hf(get_cell(nb_first_child.forward(nb_children[i])), nb_children_faces[i], hf_handle_t());
	}
}

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
void parallel_mesh<CT0, CT1, CT2, CT3>::print_cell_to_first_child(std::ostream& out) const
{
	out << "local cell-to-first-child (size = " << _local_cells.size() << ", used size = " << _local_cells.used_size() << ")" << std::endl;
	for (auto it = _local_cells.cbegin(); it != _local_cells.cend(); ++it)
	{
		cell_handle_t first_child = parallel_mesh<CT0, CT1, CT2, CT3>::first_child(*it);
		out << it.at_index() << ": ";
		if (first_child.is_valid())
		{
			if (first_child.is_ghost())	out << "g";
			out << first_child.index() << std::endl;
		}
		else
			out << "--" << std::endl;
	}
	out << "ghost cell-to-first-child (size = " << _ghost_cells.size() << ", used size = " << _ghost_cells.used_size() << ")" << std::endl;
	for (auto it = _ghost_cells.cbegin(); it != _ghost_cells.cend(); ++it)
	{
		cell_handle_t first_child = parallel_mesh<CT0, CT1, CT2, CT3>::first_child(*it);
		out << it.at_index() << ": ";
		if (first_child.is_valid())
		{
			if (first_child.is_ghost())	out << "g";
			out << first_child.index() << std::endl;
		}
		else
			out << "--" << std::endl;
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
