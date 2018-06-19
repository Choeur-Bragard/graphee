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
  /** Constructor of the Pagerank object
   * 
   * @param properties The pointer to the properties of the graph
   * @param adjency_matrix The pointer to the adjacency matrix of the graph
   * @param damp The damping factor for the pagerank (by default 0.85)
   *
   */
  Pagerank (Properties* properties, DiskSparseMatrixT* adjency_matrix, float damp = 0.85) :
    props(properties), adj_mat(adjency_matrix), damp(damp) {}

  /** Compute the Pagerank
   *
   * @param niter Number of iterations
   */
  void compute_pagerank(uint64_t niter);

  /** Save in ASCII file the top n-th sites
   *
   * @param ntops Prints the ntops-th to the file
   * @param out_filename Name of the file to print
   */
  void save_top_pagerank(uint64_t ntops, std::string out_filename) {}

private:
  const float damp; /*!< The pagerank damping factor */
  Properties* props; /*!< Pointer to the graph properties */
  DiskSparseMatrixT* adj_mat; /*!< Pointer to the ajacency matrix */

  DiskVector<Vector<float>> pagerank; /*!< Disk vector to save pagerank */
  DiskVector<Vector<float>> pagerank_itp1; /*!< Temporary disk vector to save pagerank at iteration+1 */
  DiskVector<Vector<float>> out_bounds; /*!< Disk vector to save the number of out links */
}; // class graphe::Pagerank

template <typename DiskSparseMatrixT>
void Pagerank<DiskSparseMatrixT>::compute_pagerank(uint64_t niters)
{
  /** One verifies that the ajacency matrix is not empty
   */
  if (adj_mat->empty())
  {
    print_error("Cannot compute PageRank, the \'adjency_matrix\' is empty");
    exit(-1);
  }

  /** 
   * We compute the number of out links for each page
   * One uses a trick to declare as low as possible disk vector.
   * Here, for instance we use temporary the 'pr' disk vector.
   */
  pagerank = std::move(DiskVector<Vector<float>>(props, "pr", 1.));
  out_bounds = std::move(DiskVector<Vector<float>>(props, "ob", 0.));
  out_bounds.dmat_prod_dvec(1., *adj_mat, pagerank);

  /**
   * After computing the out_bounds with the pagerank disk vector,
   * we initiate it to 1/N.
   */
  pagerank = std::move(DiskVector<Vector<float>>(props, "pr",  1./((float)props->nvertices)));
  for (uint64_t loop_id = 0; loop_id < niters; loop_id++)
  {
    /**
     * Printing the logs in strong fond (green color)
     */
    std::ostringstream oss;
    oss << "Start Pagerank loop #" << loop_id;
    print_strong_log(oss.str());

    /**
     * One inits the pagerank vector at iteration T+1:
     * \f$ PR_{t+1} = \frac{1-d}{N} \f$
     */
    pagerank_itp1 = std::move(DiskVector<Vector<float>>(props, "prp1", (1. - damp)/((float)props->nvertices)));

    /**
     * We compute \f$ PR_{t+1} += d \frac{A PR_{t}}{O} \f$
     */
    pagerank_itp1.dmat_prod_dvec_over_dvec(damp, *adj_mat, pagerank, out_bounds);

    /**
     * We swap the values between disk vector \f$ PR \f$ and \f$ PR_{t+1} \f$
     */
    pagerank_itp1.swap(pagerank);
  }
}

} // namespace graphee

#endif // GRAPEE_PAGERANK_H__
