#ifndef GRAPHEE_PAGERANK_H__
#define GRAPHEE_PAGERANK_H__

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cstring>
#include <vector>

#include "graphee.hpp"

namespace graphee
{

template <typename diskSpMatT>
class pagerank
{
public:
  pagerank (properties& properties, diskSpMatT& adjency_matrix) : props(properties), adj_mat(adjency_matrix) {}

  void compute_pagerank(uint64_t niter);

private:
  const float damping {0.85};
  properties& props;
  diskSpMatT& adj_mat;

  diskVector<vector<float>> pagerank;
  diskVector<vector<float>> pagerank_itp1;
  diskVector<vector<float>> out_bounds;
}; // class pagerank

template <typename gpe_dmat_t, typename gpe_dvec_t>
void gpePageRank<gpe_dmat_t, gpe_dvec_t>::calc_page_rank(uint64_t niter)
{
  if (adj_mat.empty())
  {
    gpe_error("Cannot compute PageRank, the \'adjency_matrix\' is empty");
    exit(-1);
  }

  // We use it as a first guess
  pagerank.set_init_value(1.);
  out_bounds.set_init_value(0.);
  out_bounds.alpha_mat_vec_prod(1., adj_mat, pagerank);

  // Necesseray for the rest of the calculation
  pagerank.set_init_value(1. / ((float)props.nslices));

  for (uint64_t loop_id = 0; loop_id < niter; loop_id++)
  {
    log.str("");
    log << "Start PageRank loop #" << loop_id;
    gpe_step(log.str());

    pagerank_itp1.set_init_value(0.);
    pagerank_itp1.alpha_mat_vec_prod(damp, adj_mat, pagerank);

    pagerank_itp1 += (1. - damp) / ((float)props.nvertices);

    pagerank_itp1.swap(pagerank);
  }
}

} // namespace graphee

#endif // GRAPEE_PAGERANK_H__
