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

#include <mpi.h>
#include <cstdlib>
#include <iostream>
#include <fstream>

#include "tests.h"

int main (int argc, char* argv[])
{
  MPI_Init(&argc, &argv);

  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  int rank_to_display = 0;
  if (argc > 1) rank_to_display = std::atoi(argv[1]);

  if (rank == 0)
  {
    if (test_mesh_construction())
      std::cout << "mesh construciton tests failed!" << std::endl << std::endl;
    else
      std::cout << "mesh construciton tests passed!" << std::endl << std::endl;
  }

  int ret = test_mesh_partitioning(rank_to_display);
  if (rank == rank_to_display)
  {
    if (ret)
      std::cout << "mesh partitioning tests failed!" << std::endl;
    else
      std::cout << "mesh partitioning tests passed!" << std::endl;
  }

  MPI_Finalize();

  return 0;
}

