#include <iostream>
#include <fstream>
#include <vector>
#include <string>

#include "gpe_props.h"
#include "gpe_diskmat.h"
#include "gpe_bsmat_csr.h"

using namespace graphee;

int main (int argc, char** argv) {
  std::string name ("cc18q1");
  gpe_props props (name, 98, 4,
      8, 15, gpe_props::GB, 256, gpe_props::MB);

  std::ifstream filelist(argv[1]);
  std::vector<std::string> filenames;
  std::string filename;
  while (filelist.good()) {
    filelist >> filename;
    filenames.push_back(filename);
  }

  gpe_diskmat<gpe_bsmat_csr<uint32_t>> dmat (props, "adj");

  dmat.load_edgelist (filenames, gpe_utils::GZ, gpe_utils::UO | gpe_utils::TRANS);
}
