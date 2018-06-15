#ifndef GRAPHEE_PAGERANK_HPP__
#define GRAPHEE_PAGERANK_HPP__

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cstring>
#include <vector>

#include "graphee.hpp"

namespace graphee
{

template <typename DiskSparseMatrixT>
class Pagerank
{
public:
  Pagerank (Properties* properties, DiskSparseMatrixT* adjency_matrix, float damp) :
    props(properties), adj_mat(adjency_matrix), damp(damp) {}
  void compute_pagerank(uint64_t niter);

private:
  const float damp;
  Properties* props;
  DiskSparseMatrixT* adj_mat;

  DiskVector<Vector<float>> pagerank;
  DiskVector<Vector<float>> pagerank_itp1;
  DiskVector<Vector<float>> out_bounds;
}; // class graphe::Pagerank

template <typename DiskSparseMatrixT>
void Pagerank<DiskSparseMatrixT>::compute_pagerank(uint64_t niters)
{
  if (adj_mat->empty())
  {
    print_error("Cannot compute PageRank, the \'adjency_matrix\' is empty");
    exit(-1);
  }

  pagerank = std::move(DiskVector<Vector<float>>(props, "pr", 1.));
  out_bounds = std::move(DiskVector<Vector<float>>(props, "ob", 0.));
  out_bounds.add_xmatvec_prod(1., *adj_mat, pagerank);

  for (uint64_t loop_id = 0; loop_id < niters; loop_id++)
  {
    std::ostringstream oss;
    oss << "Start Pagerank loop #" << loop_id;
    print_strong_log(oss.str());

    pagerank_itp1 = std::move(DiskVector<Vector<float>>(props, "prp1", 0.));

    pagerank_itp1.add_xmatvec_prod(damp, *adj_mat, pagerank);

    pagerank_itp1 += (1. - damp)/((float)props->nvertices);

    pagerank_itp1.swap(pagerank);
  }
}

} // namespace graphee

#endif // GRAPEE_PAGERANK_H__
