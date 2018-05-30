#include <iostream>
#include <fstream>
#include <vector>
#include <string>

#include "graphee.h"

using namespace graphee;

int main (int argc, char** argv) {
  std::string name ("cc18q1");
  gpe_props props (name, 98, 4, 8, 15, gpe_props::GB, 256, gpe_props::MB);

  gpe_diskmat<gpe_bsmat_csr<uint32_t>> adj (props, "adj");

  gpePageRank<gpe_diskmat<gpe_bsmat_csr<uint32_t>>, gpe_diskvec<gpe_vec<float>>> gPR (props, adj);

  gPR.calc_page_rank (10);

  return 0;

