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

#ifndef MPI_DATATYPE_TRAITS_H
#define MPI_DATATYPE_TRAITS_H 

#include <mpi.h>
#include <cstdint>

namespace pmh {

  template<typename T>
  struct mpi_datatype_traits;

  template<typename T>
  constexpr MPI_Datatype mpi_datatype = mpi_datatype_traits<T>::type;

  // specializations

  template<>
  struct mpi_datatype_traits<int>
  {
    static constexpr MPI_Datatype type = MPI_INT;
  };

  template<>
  struct mpi_datatype_traits<std::int64_t>
  {
    static constexpr MPI_Datatype type = MPI_INT64_T;
  };

  template<>
  struct mpi_datatype_traits<float>
  {
    static constexpr MPI_Datatype type = MPI_FLOAT;
  };

  template<>
  struct mpi_datatype_traits<double>
  {
    static constexpr MPI_Datatype type = MPI_DOUBLE;
  };

  // more types here ...

} // namespace pmh

#endif
