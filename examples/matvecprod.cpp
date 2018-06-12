#include <iostream>
#include <fstream>
#include <vector>
#include <string>

#include "graphee.hpp"

using namespace graphee;

int main (int argc, char** argv) {
  graphee::Properties
    props (std::string("cc18q1"), 98, 4, 8, 
    10*graphee::Properties::GB, 
    256*graphee::Properties::MB);

  graphee::DiskSparseMatrix<graphee::SparseBMatrixCSR> 
    adjacency_matrix (&props, "adj");

  float init_val = 1./((float) props.nvertices);
  graphee::DiskVector<graphee::Vector<float>>
    pagerank (&props, "pr", init_val);

  graphee::DiskVector<graphee::Vector<float>>
    pagerank_itp1 (&props, "prp1", 0.);

  pagerank_itp1.add_xmatvec_prod (1., adjacency_matrix, pagerank);

  pagerank_itp1.swap(pagerank);

  return 0;
}
