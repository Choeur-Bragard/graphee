#ifndef GPE_PAGERANK_H
#define GPE_PAGERANK_H

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cstring>
#include <vector>

#include "graphee.h"

namespace graphee {

template <typename gpe_dmat_t, typename gpe_dvec_t>
class gpePageRank {
public:
  gpePageRank (gpe_props props) : props(props) {}

  gpePageRank (gpe_props props, std::vector<std::string> filenames) :
    props(props) {
    adj_mat.load_edgelist (filenames, gpe_utils::GZ, gpe_utils::UO | gpe_utils::TRANS);
  }

  gpePageRank (gpe_props props, gpe_dmat_t& adjency_matrix) :
    props(props), adj_mat(adjency_matrix);

private:
  gpe_props props;

  gpe_dmat_t ajd_mat (props, "adj");

  gpe_dvec_t pagerank (props, "pr");
  gpe_dvec_t pagerank_itp1 (props, "prp1");
  gpe_dvec_t in_bounds (props, "ib");
  gpe_dvec_t out_bounds (props "ob");
}; // class gpePageRank

} // namespace graphee

#endif // GPE_PAGERANK_H
