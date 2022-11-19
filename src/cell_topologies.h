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

#ifndef CELL_TOPOLOGIES_H
#define CELL_TOPOLOGIES_H

#include <utility>
#include <tuple>
#include <cstring>
#include <variant>

namespace pmh {

  // the shape enum
  enum class shape_3d
  {
    HEXAHEDRON = 0,
    WEDGE,
    PYRAMID,
    TETRAHEDRON
  };

  // the following topologies all follow the CGNS convention

  class hex_topology
  {
  public:
    static shape_3d shape() { return shape_3d::HEXAHEDRON; }

    // downward incidence
    static int num_vertices() { return 8; }

    static int num_edges() { return 12; }
    // NOTE: the returned pair implicitly defines the orientation of the edge: from pair.key to pair.value
    static std::pair<int, int> vertices_of_edge(int edge)
      { return std::make_pair(_edge_to_vertices[edge][0], _edge_to_vertices[edge][1]); }

    static int num_faces() { return 6; }
    static int num_vertices_of_face(int face) { return 4; }
    // NOTE: the returned vertices must be ordered counter-clockwise w.r.t. to the face's outward normal vector
    static void vertices_of_face(int face, int* out) { std::memcpy(out, _face_to_vertices[face], 4 * sizeof(int)); }

    static int num_edges_of_face(int face) { return 4; }
    // NOTE: the returned edges must be ordered counter-clockwise w.r.t. to the face's outward normal vector
    // NOTE: and are consistent with the order of vertices_of_face() above
    static void edges_of_face(int face, int* out) { std::memcpy(out, _face_to_edges[face], 4 * sizeof(int)); }

    // upward incidence
    static int num_incident_edges(int vertex) { return 3; };
    static void incident_edges(int vertex, int* out) { std::memcpy(out, _vertex_to_edges[vertex], 3 * sizeof(int)); }

    static int num_incident_faces(int vertex) { return 3; };
    static void incident_faces(int vertex, int* out) { std::memcpy(out, _vertex_to_faces[vertex], 3 * sizeof(int)); }

    // NOTE: the returned pair implicitly defines the orientation of the faces with respect to that of the edge:
    // NOTE: pair.key is the face whose outward normal is consistent with the edge's orientation defined in
    // NOTE: vertices_of_edge() by the right-hand rule; pair.value is the face whose outward normal is opposite
    // NOTE: to the edge's orientation defined in vertices_of_edge() by the right-hand rule
    static std::pair<int, int> incident_faces(int edge)
      { return std::make_pair(_edge_to_faces[edge][0], _edge_to_faces[edge][1]); }

private:
    static constexpr int _face_to_edges[6][4] { { 3, 2, 1, 0 },
                                                { 0, 5, 8, 4 },
                                                { 1, 6, 9, 5 },
                                                { 2, 7, 10, 6 },
                                                { 4, 11, 7, 3 },
                                                { 8, 9, 10, 11 } };

    static constexpr int _face_to_vertices[6][4] { { 0, 3, 2, 1 },
                                                   { 0, 1, 5, 4 },
                                                   { 1, 2, 6, 5 },
                                                   { 2, 3, 7, 6 },
                                                   { 0, 4, 7, 3 },
                                                   { 4, 5, 6, 7 } };

    static constexpr int _edge_to_vertices[12][2] { {0, 1},
                                                    {1, 2},
                                                    {2, 3},
                                                    {3, 0},
                                                    {0, 4},
                                                    {1, 5},
                                                    {2, 6},
                                                    {3, 7},
                                                    {4, 5},
                                                    {5, 6},
                                                    {6, 7},
                                                    {7, 4} };

    static constexpr int _vertex_to_edges[8][3] { {0, 4, 3},
                                                  {0, 1, 5},
                                                  {1, 2, 6},
                                                  {2, 3, 7},
                                                  {4, 8, 11},
                                                  {5, 9, 8},
                                                  {6, 10, 9},
                                                  {7, 11, 10} };

    static constexpr int _vertex_to_faces[8][3] { {0, 1, 4},
                                                  {0, 2, 1},
                                                  {0, 3, 2},
                                                  {0, 4, 3},
                                                  {1, 5, 4},
                                                  {2, 5, 1},
                                                  {2, 3, 5},
                                                  {3, 4, 5} };

    static constexpr int _edge_to_faces[12][2] { {1, 0},
                                                 {2, 0},
                                                 {3, 0},
                                                 {4, 0},
                                                 {4, 1},
                                                 {1, 2},
                                                 {2, 3},
                                                 {3, 4},
                                                 {5, 1},
                                                 {5, 2},
                                                 {5, 3},
                                                 {5, 4} };
  };

  class wdg_topology
  {
  public:
    static shape_3d shape() { return shape_3d::WEDGE; }

    // downward incidence
    static int num_vertices() { return 6; }

    static int num_edges() { return 9; }
    // NOTE: the returned pair implicitly defines the orientation of the edge: from pair.key to pair.value
    static std::pair<int, int> vertices_of_edge(int edge)
      { return std::make_pair(_edge_to_vertices[edge][0], _edge_to_vertices[edge][1]); }

    static int num_faces() { return 5; }
    static int num_vertices_of_face(int face) { return face < 3 ? 4 : 3; }
    // NOTE: the returned vertices must be ordered counter-clockwise w.r.t. to the face's outward normal vector
    static void vertices_of_face(int face, int* out) { std::memcpy(out, _face_to_vertices[face], 4 * sizeof(int)); }

    static int num_edges_of_face(int face) { return face < 3 ? 4 : 3; }
    // NOTE: the returned edges must be ordered counter-clockwise w.r.t. to the face's outward normal vector
    // NOTE: and are consistent with the order of vertices_of_face() above
    static void edges_of_face(int face, int* out) { std::memcpy(out, _face_to_edges[face], 4 * sizeof(int)); }

    // upward incidence
    static int num_incident_edges(int vertex) { return 3; };
    static void incident_edges(int vertex, int* out) { std::memcpy(out, _vertex_to_edges[vertex], 3 * sizeof(int)); }

    static int num_incident_faces(int vertex) { return 3; };
    static void incident_faces(int vertex, int* out) { std::memcpy(out, _vertex_to_faces[vertex], 3 * sizeof(int)); }

    // NOTE: the returned pair implicitly defines the orientation of the faces with respect to that of the edge:
    // NOTE: pair.key is the face whose outward normal is consistent with the edge's orientation defined in
    // NOTE: vertices_of_edge() by the right-hand rule; pair.value is the face whose outward normal is opposite
    // NOTE: to the edge's orientation defined in vertices_of_edge() by the right-hand rule
    static std::pair<int, int> incident_faces(int edge)
      { return std::make_pair(_edge_to_faces[edge][0], _edge_to_faces[edge][1]); }

  private:
    static constexpr int _face_to_edges[5][4] { { 0, 4, 6, 3 },
                                                { 1, 5, 7, 4 },
                                                { 2, 3, 8, 5 },
                                                { 2, 1, 0, -1 },
                                                { 6, 7, 8, -1 } };

    static constexpr int _face_to_vertices[5][4] { { 0, 1, 4, 3 },
                                                   { 1, 2, 5, 4 },
                                                   { 2, 0, 3, 5 },
                                                   { 0, 2, 1, -1 },
                                                   { 3, 4, 5, -1 } };

    static constexpr int _edge_to_vertices[9][2] { {0, 1},
                                                   {1, 2},
                                                   {2, 0},
                                                   {0, 3},
                                                   {1, 4},
                                                   {2, 5},
                                                   {3, 4},
                                                   {4, 5},
                                                   {5, 3} };

    static constexpr int _vertex_to_edges[6][3] { {0, 3, 2},
                                                  {1, 4, 0},
                                                  {2, 5, 1},
                                                  {6, 8, 3},
                                                  {7, 6, 4},
                                                  {8, 7, 5} };

    static constexpr int _vertex_to_faces[6][3] { {3, 0, 2},
                                                  {3, 1, 0},
                                                  {3, 2, 1},
                                                  {4, 2, 0},
                                                  {4, 0, 1},
                                                  {4, 1, 2} };

    static constexpr int _edge_to_faces[9][2] { {0, 3},
                                                {1, 3},
                                                {2, 3},
                                                {2, 0},
                                                {0, 1},
                                                {1, 2},
                                                {4, 0},
                                                {4, 1},
                                                {4, 2} };
  };

  class prm_topology
  {
  public:
    static shape_3d shape() { return shape_3d::PYRAMID; }

    // downward incidence
    static int num_vertices() { return 5; }

    static int num_edges() { return 8; }
    // NOTE: the returned pair implicitly defines the orientation of the edge: from pair.key to pair.value
    static std::pair<int, int> vertices_of_edge(int edge)
      { return std::make_pair(_edge_to_vertices[edge][0], _edge_to_vertices[edge][1]); }

    static int num_faces() { return 5; }
    static int num_vertices_of_face(int face) { return face == 0 ? 4 : 3; }
    // NOTE: the returned vertices must be ordered counter-clockwise w.r.t. to the face's outward normal vector
    static void vertices_of_face(int face, int* out) { std::memcpy(out, _face_to_vertices[face], 4 * sizeof(int)); }

    static int num_edges_of_face(int face) { return face == 0 ? 4 : 3; }
    // NOTE: the returned edges must be ordered counter-clockwise w.r.t. to the face's outward normal vector
    // NOTE: and are consistent with the order of vertices_of_face() above
    static void edges_of_face(int face, int* out) { std::memcpy(out, _face_to_edges[face], 4 * sizeof(int)); }

    // upward incidence
    static int num_incident_edges(int vertex) { return vertex < 4 ? 3 : 4; };
    static void incident_edges(int vertex, int* out) { std::memcpy(out, _vertex_to_edges[vertex], 4 * sizeof(int)); }

    static int num_incident_faces(int vertex) { return vertex < 4 ? 3 : 4; };
    static void incident_faces(int vertex, int* out) { std::memcpy(out, _vertex_to_faces[vertex], 4 * sizeof(int)); }

    // NOTE: the returned pair implicitly defines the orientation of the faces with respect to that of the edge:
    // NOTE: pair.key is the face whose outward normal is consistent with the edge's orientation defined in
    // NOTE: vertices_of_edge() by the right-hand rule; pair.value is the face whose outward normal is opposite
    // NOTE: to the edge's orientation defined in vertices_of_edge() by the right-hand rule
    static std::pair<int, int> incident_faces(int edge)
      { return std::make_pair(_edge_to_faces[edge][0], _edge_to_faces[edge][1]); }

private:
    static constexpr int _face_to_edges[5][4] { { 3, 2, 1, 0 },
                                                { 0, 5, 4, -1 },
                                                { 1, 6, 5, -1 },
                                                { 2, 7, 6, -1 },
                                                { 3, 4, 7, -1 } };

    static constexpr int _face_to_vertices[5][4] { { 0, 3, 2, 1 },
                                                   { 0, 1, 4, -1 },
                                                   { 1, 2, 4, -1 },
                                                   { 2, 3, 4, -1 },
                                                   { 3, 0, 4, -1 } };

    static constexpr int _edge_to_vertices[8][2] { {0, 1},
                                                   {1, 2},
                                                   {2, 3},
                                                   {3, 0},
                                                   {0, 4},
                                                   {1, 4},
                                                   {2, 4},
                                                   {3, 4} };

    static constexpr int _vertex_to_edges[5][4] { {0, 4, 3, -1},
                                                  {1, 5, 0, -1},
                                                  {2, 6, 1, -1},
                                                  {3, 7, 2, -1},
                                                  {4, 5, 6, 7} };

    static constexpr int _vertex_to_faces[5][4] { {0, 1, 4, -1},
                                                  {0, 2, 1, -1},
                                                  {0, 3, 2, -1},
                                                  {0, 4, 3, -1},
                                                  {1, 2, 3, 4} };

    static constexpr int _edge_to_faces[8][2] { {1, 0},
                                                {2, 0},
                                                {3, 0},
                                                {4, 0},
                                                {4, 1},
                                                {1, 2},
                                                {2, 3},
                                                {3, 4} };
  };

  class tet_topology
  {
  public:
    static shape_3d shape() { return shape_3d::TETRAHEDRON; }

    // downward incidence
    static int num_vertices() { return 4; }

    static int num_edges() { return 6; }
    // NOTE: the returned pair implicitly defines the orientation of the edge: from pair.key to pair.value
    static std::pair<int, int> vertices_of_edge(int edge)
      { return std::make_pair(_edge_to_vertices[edge][0], _edge_to_vertices[edge][1]); }

    static int num_faces() { return 4; }
    static int num_vertices_of_face(int face) { return 3; }
    // NOTE: the returned vertices must be ordered counter-clockwise w.r.t. to the face's outward normal vector
    static void vertices_of_face(int face, int* out) { std::memcpy(out, _face_to_vertices[face], 3 * sizeof(int)); }

    static int num_edges_of_face(int face) { return 3; }
    // NOTE: the returned edges must be ordered counter-clockwise w.r.t. to the face's outward normal vector
    // NOTE: and are consistent with the order of vertices_of_face() above
    static void edges_of_face(int face, int* out) { std::memcpy(out, _face_to_edges[face], 3 * sizeof(int)); }

    // upward incidence
    static int num_incident_edges(int vertex) { return 3; };
    static void incident_edges(int vertex, int* out) { std::memcpy(out, _vertex_to_edges[vertex], 3 * sizeof(int)); }

    static int num_incident_faces(int vertex) { return 3; };
    static void incident_faces(int vertex, int* out) { std::memcpy(out, _vertex_to_faces[vertex], 3 * sizeof(int)); }

    // NOTE: the returned pair implicitly defines the orientation of the faces with respect to that of the edge:
    // NOTE: pair.key is the face whose outward normal is consistent with the edge's orientation defined in
    // NOTE: vertices_of_edge() by the right-hand rule; pair.value is the face whose outward normal is opposite
    // NOTE: to the edge's orientation defined in vertices_of_edge() by the right-hand rule
    static std::pair<int, int> incident_faces(int edge)
      { return std::make_pair(_edge_to_faces[edge][0], _edge_to_faces[edge][1]); }

  private:
    static constexpr int _face_to_edges[4][3] { { 2, 1, 0 },
                                                { 0, 4, 3 },
                                                { 1, 5, 4 },
                                                { 2, 3, 5 } };

    static constexpr int _face_to_vertices[4][3] { { 0, 2, 1 },
                                                   { 0, 1, 3 },
                                                   { 1, 2, 3 },
                                                   { 2, 0, 3 } };

    static constexpr int _edge_to_vertices[6][2] { {0, 1},
                                                   {1, 2},
                                                   {2, 0},
                                                   {0, 3},
                                                   {1, 3},
                                                   {2, 3} };

    static constexpr int _vertex_to_edges[4][3] { {0, 2, 3},
                                                  {0, 1, 4},
                                                  {1, 2, 5},
                                                  {3, 4, 5} };

    static constexpr int _vertex_to_faces[4][3] { {0, 1, 3},
                                                  {0, 1, 2},
                                                  {0, 2, 3},
                                                  {1, 2, 3} };

    static constexpr int _edge_to_faces[6][2] { {1, 0},
                                                {2, 0},
                                                {3, 0},
                                                {3, 1},
                                                {1, 2},
                                                {2, 3} };
  };


  // extract topology type from cell type
  template<typename T>
  using topology_t = typename T::topology;

  // CT - single cell type
  template<typename CT>
  class single_shape_topology
  {
  public:
    static shape_3d shape(const CT& cell)
      { return topo_type::shape(); }

    static int num_vertices(const CT& cell)
      { return topo_type::num_vertices(); }

    static int num_edges(const CT& cell)
      { return topo_type::num_edges(); }

    static std::pair<int, int> vertices_of_edge(const CT& cell, int edge)
      { return topo_type::vertices_of_edge(edge); }

    static int num_faces(const CT& cell)
      { return topo_type::num_faces(); }

    static int num_vertices_of_face(const CT& cell, int face)
      { return topo_type::num_vertices_of_face(face); }

    static void vertices_of_face(const CT& cell, int face, int* out)
      { return topo_type::vertices_of_face(face, out); }

    static int num_edges_of_face(const CT& cell, int face)
      { return topo_type::num_edges_of_face(face); }

    static void edges_of_face(const CT& cell, int face, int* out)
      { return topo_type::edges_of_face(face, out); }

    static int num_incident_edges(const CT& cell, int vertex)
      { return topo_type::num_incident_edges(vertex); }

    static void incident_edges(const CT& cell, int vertex, int* out)
      { return topo_type::incident_edges(vertex, out); }

    static int num_incident_faces(const CT& cell, int vertex)
      { return topo_type::num_incident_faces(vertex); }

    static void incident_faces(const CT& cell, int vertex, int* out)
      { return topo_type::incident_faces(vertex, out); }
    
    // NOTE: the returned pair implicitly defines the orientation of the faces with respect to that of the edge:
    // NOTE: pair.key is the face whose outward normal is consistent with the edge's orientation defined in
    // NOTE: vertices_of_edge() by the right-hand rule; pair.value is the face whose outward normal is opposite
    // NOTE: to the edge's orientation defined in vertices_of_edge() by the right-hand rule
    static std::pair<int, int> incident_faces(const CT& cell, int edge)
      { return topo_type::incident_faces(edge); }

  private:
    using topo_type = topology_t<CT>;
  };

  // CTS - multiple cell types
  template<typename... CTS>
  class mixed_shape_topology
  {
  public:
    using cell_type = std::variant<CTS...>;

    static shape_3d shape(const cell_type& cell)
      { return std::visit([](auto&& t) { return t.shape(); }, _topos[cell.index()]); }

    static int num_vertices(const cell_type& cell)
      { return std::visit([](auto&& t) { return t.num_vertices(); }, _topos[cell.index()]); }

    static int num_edges(const cell_type& cell)
      { return std::visit([](auto&& t) { return t.num_edges(); }, _topos[cell.index()]); }

    static std::pair<int, int> vertices_of_edge(const cell_type& cell, int edge)
      { return std::visit([edge](auto&& t) { return t.vertices_of_edge(edge); }, _topos[cell.index()]); }

    static int num_faces(const cell_type& cell)
      { return std::visit([](auto&& t) { return t.num_faces(); }, _topos[cell.index()]); }

    static int num_vertices_of_face(const cell_type& cell, int face)
      { return std::visit([face](auto&& t) { return t.num_vertices_of_face(face); }, _topos[cell.index()]); }

    static void vertices_of_face(const cell_type& cell, int face, int* out)
      { std::visit([face, out](auto&& t) { return t.vertices_of_face(face, out); }, _topos[cell.index()]); }

    static int num_edges_of_face(const cell_type& cell, int face)
      { return std::visit([face](auto&& t) { return t.num_edges_of_face(face); }, _topos[cell.index()]); }

    static void edges_of_face(const cell_type& cell, int face, int* out)
      { std::visit([face, out](auto&& t) { return t.edges_of_face(face, out); }, _topos[cell.index()]); }

    static int num_incident_edges(const cell_type& cell, int vertex)
      { return std::visit([vertex](auto&& t) { return t.num_incident_edges(vertex); }, _topos[cell.index()]); }

    static void incident_edges(const cell_type& cell, int vertex, int* out)
      { std::visit([vertex, out](auto&& t) { return t.incident_edges(vertex, out); }, _topos[cell.index()]); }

    static int num_incident_faces(const cell_type& cell, int vertex)
      { return  std::visit([vertex](auto&& t) { return t.num_incident_faces(vertex); }, _topos[cell.index()]); }

    static void incident_faces(const cell_type& cell, int vertex, int* out)
      { std::visit([vertex, out](auto&& t) { return t.incident_faces(vertex, out); }, _topos[cell.index()]); }
    
    // NOTE: the returned pair implicitly defines the orientation of the faces with respect to that of the edge:
    // NOTE: pair.key is the face whose outward normal is consistent with the edge's orientation defined in
    // NOTE: vertices_of_edge() by the right-hand rule; pair.value is the face whose outward normal is opposite
    // NOTE: to the edge's orientation defined in vertices_of_edge() by the right-hand rule
    static std::pair<int, int> incident_faces(const cell_type& cell, int edge)
      { return std::visit([edge](auto&& t) { return t.incident_faces(edge); }, _topos[cell.index()]); }

  private:
    // construct variant of topologies corresponding to the given cell types
    using topo_type = std::variant<topology_t<CTS>...>;
    static const topo_type _topos[sizeof...(CTS)];
  };

} // namespace pmh

#endif
