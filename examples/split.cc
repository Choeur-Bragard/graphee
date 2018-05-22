#include <iostream>
#include <fstream>
#include <vector>
#include <string>

#include "gpe_props.h"
#include "gpe_diskmat.h"
#include "gpe_bsmat_csr.h"

using namespace graphee;

int main (int argc, char** argv) {
  std::string name ("large");
  gpe_props props (name, 3564, 10,
      80, 400, gpe_props::GB, 1, gpe_props::GB);

  std::ifstream filelist(argv[1]);
  std::vector<std::string> filenames;
  std::string filename;
  while (filelist.good()) {
    filelist >> filename;
    filenames.push_back(filename);
  }

  gpe_diskmat<gpe_bsmat_csr<uint32_t>, uint32_t> dmat (props);

  dmat.load_edgelist (filenames, gpe_utils::GZ, gpe_utils::UO | gpe_utils::TRANS);
}
