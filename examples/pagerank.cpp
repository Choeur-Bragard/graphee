#include <iostream>
#include <fstream>
#include <vector>
#include <string>

#include "graphee.hpp"

using namespace graphee;

int main (int argc, char** argv)
{
  graphee::Properties
  props (std::string("cc18q1"), 98, 4, 8,
         10*graphee::Properties::GB,
         256*graphee::Properties::MB);

  graphee::DiskSparseMatrix<graphee::SparseBMatrixCSR>
  adjacency_matrix (&props, "adj");

  graphee::Pagerank<graphee::DiskSparseMatrix<graphee::SparseBMatrixCSR>>
    pagerank(&props, &adjacency_matrix, 0.85);

  pagerank.compute_pagerank(10);

  return 0;
}
