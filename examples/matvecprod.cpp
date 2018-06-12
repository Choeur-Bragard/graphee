#include <iostream>
#include <fstream>
#include <vector>
#include <string>

#include "graphee.hpp"

using namespace graphee;

int main (int argc, char** argv) {
  graphee::properties
    props (std::string("cc18q1"), 98, 4, 8, 
    10*graphee::properties::GB, 
    256*graphee::properties::MB);

  graphee::diskSparseMatrix<graphee::sparseBMatrixCSR> 
    adjency_matrix (&props, "adj");

  float init_val = 1./((float) props.nvertices);
  graphee::diskVector<graphee::vector<float>>
    pagerank (&props, "pr", init_val);

  graphee::diskVector<graphee::vector<float>>
    pagerank_itp1 (&props, "prp1", 0.);

  pagerank_itp1.add_xmatvec_prod (1., adjency_matrix, pagerank);

  pagerank_itp1.swap(pagerank);

  return 0;
}
