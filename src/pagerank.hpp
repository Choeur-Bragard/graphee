#ifndef GRAPHEE_PAGERANK_HPP__
#define GRAPHEE_PAGERANK_HPP__

#include <cstring>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include "graphee.hpp"

namespace graphee {

template <typename DiskSparseMatrixT> class Pagerank {
public:
  /** Constructor of the Pagerank object
   *
   * Please take a look to: https://en.wikipedia.org/wiki/PageRank#Iterative
   * for a more precise description of the pagerank algorithm.
   *
   * @param properties The pointer to the properties of the graph
   * @param adjency_matrix The pointer to the adjacency matrix of the graph
   * @param damp The damping factor for the pagerank (by default 0.85)
   *
   */
  Pagerank(Properties *properties, DiskSparseMatrixT *adjency_matrix,
           double damp = 0.85)
      : props(properties), adj_mat(adjency_matrix), damp(damp) {}

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
  const double damp;           ///< The pagerank damping factor
  Properties *props;          ///< Pointer to the graph properties
  DiskSparseMatrixT *adj_mat; ///< Pointer to the ajacency matrix

  DiskVector<Vector<double>> pagerank; ///< Disk vector to save pagerank
  DiskVector<Vector<double>>
      pagerank_itp1; ///< Temporary disk vector to save pagerank at iteration+1
  DiskVector<Vector<double>>
      out_bounds; ///< Disk vector to save the number of out links
  
  void pagerankStatistic(DiskVector<Vector<double>> &pagerank_itp1,
                       DiskVector<Vector<double>> &pagerank,
                       DiskVector<Vector<double>> &out_bounds, double &sinkScore,
                       double &sumScore, double &scoreVariation);
};

template <typename DiskSparseMatrixT>
void Pagerank<DiskSparseMatrixT>::compute_pagerank(uint64_t niters) {
  // One verifies that the ajacency matrix is not empty
  if (adj_mat->empty()) {
    print_error("Cannot compute PageRank, the \'adjency_matrix\' is empty");
    exit(-1);
  }

  // number of outlinks of each node
  out_bounds = std::move(DiskVector<Vector<double>>(props, "ob", 0));
  out_bounds.dmat_columns_sum(*adj_mat);

  // pagerank is initiated to 1/N.
  pagerank = std::move(
      DiskVector<Vector<double>>(props, "pr", 1. / ((double)props->nvertices)));

  // count the number of nodes without outlinks
  uint64_t n_sink_nodes = out_bounds.countZeros();
  // amount of score given by the sink nodes to each nodes
  double sinkScore = 1. / ((double)props->nvertices) * n_sink_nodes;
  // pagerank variation with the previous iteration (arbirtraily initialised at
  // 0)
  double scoreVariation = 0;
  // sum of all pagerank score
  double sumScore = 1.;

  for (uint64_t loop_id = 0; loop_id < niters; loop_id++) {
    std::ostringstream oss;
    oss << "Start Pagerank loop #" << loop_id;
    print_strong_log(oss.str());

    // One inits the pagerank vector at iteration T+1 with the score received by
    // random jump plus the redistribution from sink node.
    pagerank_itp1 = std::move(DiskVector<Vector<double>>(
        props, "prp1",
        (1 - sinkScore) * (1. - damp) / ((double)props->nvertices) +
            sinkScore / ((double)props->nvertices)));

    // We compute \f$ PR_{t+1} \f$
    pagerank_itp1.dmat_prod_dvec_over_dvec(damp, *adj_mat, pagerank,
                                           out_bounds);

    // We compute statistics
    pagerankStatistic(pagerank_itp1, pagerank, out_bounds, sinkScore, sumScore,
                      scoreVariation);

    // Display statistics report
    oss.str("");
    oss.clear();
    oss << "scoreVariation : " << scoreVariation;
    print_strong_log(oss.str());
    oss.str("");
    oss.clear();
    oss << "sumScore : " << sumScore;
    print_strong_log(oss.str());
    oss.str("");
    oss.clear();
    oss << "sinkScore : " << sinkScore;
    print_strong_log(oss.str());
    oss.str("");
    oss.clear();
    oss << "end of iteration #" << loop_id;
    print_strong_log(oss.str());

    pagerank_itp1.swap(pagerank);
  }
}

/** Compute Pagerank statistics
 *
 * @param scoreVariation : variation of score with the previous iteration. can
 * be used to get an estimation of the convergence or as a stopping criterion.
 * @param sinSckore : amount of score held by nodes without outlinks. this
 * statistic is required for the next pagerank iteration.
 * @param sumScore : sum of each node's score. value far from
 * one indicate a bug or numericale error.
 */
template <typename DiskSparseMatrixT>
void Pagerank<DiskSparseMatrixT>::pagerankStatistic(DiskVector<Vector<double>> &pagerank_itp1,
                       DiskVector<Vector<double>> &pagerank,
                       DiskVector<Vector<double>> &out_bounds, double &sinkScore,
                       double &sumScore, double &scoreVariation) {

  sinkScore = 0;
  sumScore = 0;
  scoreVariation = 0;

#pragma omp parallel for reduction(+ : sinkScore, sumScore, scoreVariation)
  for (uint64_t slice = 0; slice < pagerank.get_nslices(); slice++) {
    Vector<double> out_bounds_vec(std::move(out_bounds.get_slice(slice)));
    Vector<double> pagerank_vec(std::move(pagerank.get_slice(slice)));
    Vector<double> pagerank_itp1_vec(std::move(pagerank_itp1.get_slice(slice)));

    for (int i = 0; i < pagerank_vec.get_lines(); i++) {
      sumScore += pagerank_itp1_vec[i];
      scoreVariation += (pagerank_vec[i] - pagerank_itp1_vec[i]) *
                        (pagerank_vec[i] - pagerank_itp1_vec[i]);

      if (out_bounds_vec[i] == 0) {
        sinkScore += pagerank_itp1_vec[i];
      }
    }
  }
}

} // namespace graphee

#endif // GRAPEE_PAGERANK_H__
