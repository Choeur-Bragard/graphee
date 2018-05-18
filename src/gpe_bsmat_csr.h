#ifndef GPE_BSMAT_CSR_H
#define GPE_BSMAT_CSR_H

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

#include "gpe_props.h"
#include "gpe_utils.h"

namespace graphee {

/*! \brief Boolean sparse matrix
 *
 * It defines a sparse matrix filled only with elements,
 * of values \{0,1\}. This is the case of a non-weighted-edges
 * graph.
 */

template <class idx_t>
class gpe_bsmat_csr {
public:
  enum {BIN, SNAPPY};

  gpe_bsmat_csr (gpe_props* i_prop);
  gpe_bsmat_csr (gpe_props* i_prop, idx_t i_m, idx_t i_nnz);
  ~gpe_bsmat_csr ();

  void ordered_add (idx_t i, idx_t j);

  void insert (idx_t i, idx_t j);
  void remove (idx_t i, idx_t j);

  void save (std::string name, int fileformat = BIN, int64_t offl = 0, int64_t offc = 0);
  void load (std::string name, int fileformat = BIN);

private:
  gpe_props *prop;

  bool is_alloc {false};

  idx_t m;
  idx_t nnz;

  int64_t offl;
  int64_t offc;

  idx_t ordpos {0};

  idx_t *ia;
  idx_t *ja;
}; // class gpe_bsmat_csr

template <class idx_t>
gpe_bsmat_csr<idx_t>::gpe_bsmat_csr (gpe_props* i_prop) {
  prop = i_prop;
}

template <class idx_t>
gpe_bsmat_csr<idx_t>::gpe_bsmat_csr (gpe_props* i_prop, idx_t i_m, idx_t i_nnz) {
  prop = i_prop;
  m = i_m;
  nnz = i_nnz;

  if ((nnz+m+1)*sizeof(idx_t) < prop->ram_limit) {
    ia = new idx_t [m+1];
    ja = new idx_t [nnz];
    is_alloc = true;
  } else {
    std::cerr << "[GRAPHEE] [GPE_BSMAT_CSR] Requested size is beyond RAMLIMIT" << std::end;
    exit (-1);
  }
}

template <class idx_t>
gpe_bsmat_csr<idx_t>::~gpe_bsmat_csr () {
  delete[] ia;
  delete[] ja;
}

template <class idx_t>
void gpe_bsmat_csr<idx_t>::ordered_add (idx_t i, idx_t j) {
}

} // namespace graphee

#endif // GPE_BSMAT_CSR_H
