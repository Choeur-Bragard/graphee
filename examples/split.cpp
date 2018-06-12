#include <iostream>
#include <fstream>
#include <vector>
#include <string>

#include "graphee.hpp"

int main (int argc, char** argv) {
  graphee::Properties
    props (std::string("cc18q1"), 98, 4, 8, 
    10*graphee::Properties::GB, 
    256*graphee::Properties::MB);

  std::ifstream filelist(argv[1]);
  std::vector<std::string> filenames;
  std::string filename;
  while (filelist.good()) {
    filelist >> filename;
    filenames.push_back(filename);
  }
  filenames.pop_back();
  filelist.close();

  graphee::DiskSparseMatrix<graphee::SparseBMatrixCSR> 
    adjacency_matrix (&props, "adj");

  adjacency_matrix.load_edgelist(filenames);
}
