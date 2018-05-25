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
  gpe_props props (name, 98, 4,
      8, 15, gpe_props::GB, 256, gpe_props::MB);

  gpe_diskmat<gpe_bsmat_csr<uint32_t>, uint32_t> dmat (props);

  gpe_bsmat_csr<uint32_t> smat (props);

  dmat.get_matrix_block (0, 0, smat);

  if (smat.verify()) {
    std::cout << "Correct loading of the matrix" << std::endl;
    return 0;
  } else {
    std::cout << "Incorrect loading of the matrix" << std::endl;
    return -1;
  }
}
