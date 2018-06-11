#include <iostream>
#include <fstream>
#include <vector>
#include <string>

#include "graphee.h"

int main (int argc, char** argv) {
  graphee::properties
    props (std::string("cc18q1"), 98, 4, 8, 
    10*graphee::properties::GB, 
    256*graphee::properties::MB);

  std::ifstream filelist(argv[1]);
  std::vector<std::string> filenames;
  std::string filename;
  while (filelist.good()) {
    filelist >> filename;
    filenames.push_back(filename);
  }
  filenames.pop_back();
  filelist.close();

  graphee::diskSparseMatrix<graphee::sparseBMatrixCSR> 
    adjency_matrix (&props, "adj");

  adjency_matrix.load_edgelist(filenames);
}
