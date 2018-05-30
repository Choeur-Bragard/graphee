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
    props(props), adj_mat(adjency_matrix) {};

  void calc_page_rank (uint64_t niter);

private:
  gpe_props props;

  gpe_dmat_t adj_mat;

  gpe_dvec_t pagerank;
  gpe_dvec_t pagerank_itp1;
  gpe_dvec_t out_bounds;

  float damp {0.85};
}; // class gpePageRank

template <typename gpe_dmat_t, typename gpe_dvec_t>
void gpePageRank<gpe_dmat_t, gpe_dvec_t>::calc_page_rank (uint64_t niter) {
  if (adj_mat.empty()) {
    gpe_error ("Cannot compute PageRank, the \'adjency_matrix\' is empty");
    exit (-1);
  }

  // We use it as a first guess 
  pagerank.set_init_value (1.);
  out_bounds.mat_vec_prod (adj_mat, pagerank);

  // Necesseray for the rest of the calculation
  pagerank.set_init_value (1./((float) props.nslices));

  for (uint64_t loop_id = 0; loop_id < niter; loop_id++) {
    pagerank_itp1.init_val (0.);
    pagerank_itp1.alpha_mat_vec_prod (damp, adj_mat, pagerank);

    pagerank_itp1 += (1.-damp)/((float) props.nvertices);

    pagerank_itp1.swap (pagerank);
  }
}

} // namespace graphee

#endif // GPE_PAGERANK_H
