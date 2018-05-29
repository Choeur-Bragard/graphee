#include <iostream>
#include <fstream>
#include <vector>
#include <string>

#include "gpe_props.h"
#include "gpe_diskmat.h"
#include "gpe_bsmat_csr.h"
#include "gpe_utils.h"

using namespace graphee;

int main (int argc, char** argv) {
  std::string name ("cc18q1");
  gpe_props props (name, 98, 4,
      8, 15, gpe_props::GB, 256, gpe_props::MB);

  gpe_diskmat<gpe_bsmat_csr<uint32_t>> dmat (props, "adj");

  gpe_bsmat_csr<uint32_t> smat (props);

  dmat.get_matrix_block (0, 0, smat);

  if (smat.verify()) {
    gpe_log ("Correct loading of the matrix");
    return 0;
  } else {
    gpe_error ("Incorrect loading of the matrix");
    return -1;
  }
}
