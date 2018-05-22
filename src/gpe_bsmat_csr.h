#ifndef GPE_BSMAT_CSR_H
#define GPE_BSMAT_CSR_H

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cstring>
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

  gpe_bsmat_csr (gpe_props& i_prop);
  gpe_bsmat_csr (gpe_props& i_prop, idx_t i_m, idx_t i_nnz);
  ~gpe_bsmat_csr ();

  void sorted_fill (idx_t i, idx_t j);

  void set_offsets (uint64_t offl, uint64_t offc);

  void insert (idx_t i, idx_t j);
  void remove (idx_t i, idx_t j);

  void save (std::string name, int fileformat = BIN, int64_t offl = 0, int64_t offc = 0);
  void load (std::string name);

  size_t size ();
  bool verify ();

private:
  gpe_props prop;

  bool is_alloc {false};

  idx_t last_id {0};
  idx_t m;
  idx_t nnz;

  uint64_t offl;
  uint64_t offc;

  idx_t *ia;
  idx_t *ja;
}; // class gpe_bsmat_csr

/*! Clean constructor */
template <class idx_t>
gpe_bsmat_csr<idx_t>::gpe_bsmat_csr (gpe_props& i_prop) {
  prop = i_prop;
  m = 0;
  nnz = 0;
}

/*! Constructor for predefined matrix size */
template <class idx_t>
gpe_bsmat_csr<idx_t>::gpe_bsmat_csr (gpe_props& i_prop, idx_t i_m, idx_t i_nnz) {
  prop = i_prop;
  m = i_m;
  nnz = i_nnz;

  if ((nnz+m+1)*sizeof(idx_t) < prop.ram_limit) {
    ia = new idx_t [m+1];
    ia[0] = 0;
    ja = new idx_t [nnz];
    is_alloc = true;
  } else {
    std::cerr << "[GRAPHEE] [GPE_BSMAT_CSR] Requested size is beyond RAMLIMIT" << std::endl;
    exit (-1);
  }
}

/*! General destructor of the class */
template <class idx_t>
gpe_bsmat_csr<idx_t>::~gpe_bsmat_csr () {
  delete[] ia;
  delete[] ja;
}

/*! Standard fill of the matrix */
template <class idx_t>
void gpe_bsmat_csr<idx_t>::sorted_fill (idx_t i, idx_t j) {
  for (idx_t l = last_id+1; l <= i; l++) {
    ia[l+1] = ia[l];
  }

  ja[ia[i+1]] = j;
  ia[i+1]++;

  last_id = i;
}

template <class idx_t>
void gpe_bsmat_csr<idx_t>::set_offsets (uint64_t i_offl, uint64_t i_offc) {
  offl = i_offl;
  offc = i_offc;
}

/*! Inserting element in a CSR matrix */
template <class idx_t>
void gpe_bsmat_csr<idx_t>::insert (idx_t i, idx_t j) {
}

/*! Remove element of the CSR matrix*/
template <class idx_t>
void gpe_bsmat_csr<idx_t>::remove (idx_t i, idx_t j) {
}

/*! Save the matrix to a file
 *
 * One can determine the fileformat avail BIN or
 * SNAPPY for fast-light compression.
 */
template <class idx_t>
void gpe_bsmat_csr<idx_t>::save (std::string name, int fileformat, int64_t offl, int64_t) {
  std::ofstream matfp (name, std::ios_base::binary);

  std::string matrix_type ("GPE_BSMAT_CSR");
  size_t matrix_type_size = matrix_type.size();

  /* Save explicitly matrix properties */
  matfp.write ((const char*) &matrix_type_size, sizeof(size_t));
  matfp.write ((const char*) matrix_type.c_str(), matrix_type_size);

  /* Save fileformat {BIN, SNAPPY} */
  matfp.write ((const char*) &fileformat, sizeof(int));

  /* Printing offsets of the matrix if needed */
  matfp.write ((const char*) &offl, sizeof(int64_t));
  matfp.write ((const char*) &offc, sizeof(int64_t));

  /* Matrix dimension */
  matfp.write ((const char*) &m, sizeof(idx_t));
  matfp.write ((const char*) &nnz, sizeof(idx_t));

  if (fileformat == BIN) {
    matfp.write ((const char*) ia, (m+1)*sizeof(idx_t));
    matfp.write ((const char*) ja, (nnz)*sizeof(idx_t));

  } else if (fileformat == SNAPPY) {
    bool comp_succeed;

    char* ia_snappy;
    size_t ia_snappy_size;
    comp_succeed = compress_snappy ((const char*) ia, (m+1)*sizeof(idx_t), &ia_snappy, ia_snappy_size);

    if (!comp_succeed) {
      std::cerr << "[GRAPHEE] [GPE_BSMAT_CSR] SNAPPY compression of IA failed" << std::endl;

      delete[] ia_snappy;
      matfp.close();
      exit (-1);
    }

    matfp.write ((const char*) &ia_snappy_size, sizeof(size_t));
    matfp.write ((const char*) ia_snappy, ia_snappy_size);
    delete[] ia_snappy;

    char* ja_snappy;
    size_t ja_snappy_size;
    comp_succeed = compress_snappy ((const char*) ja, (nnz)*sizeof(idx_t), &ja_snappy, ja_snappy_size);

    if (!comp_succeed) {
      std::cerr << "[GRAPHEE] [GPE_BSMAT_CSR] SNAPPY compression of JA failed" << std::endl;

      delete[] ja_snappy;
      matfp.close();
      exit (-1);
    }

    matfp.write ((const char*) &ja_snappy_size, sizeof(size_t));
    matfp.write ((const char*) ja_snappy, ja_snappy_size);
    delete[] ja_snappy;
  }

  matfp.close();
}

/*! Read matrix from a file
 *
 * It determines the fileformat from the file itself
 */
template <class idx_t>
void gpe_bsmat_csr<idx_t>::load (std::string name) {
  std::ifstream matfp (name, std::ios_base::binary);

  /* Save explicitly matrix properties */
  size_t matrix_type_size;
  matfp.read ((char*) &matrix_type_size, sizeof(size_t));

  char matrix_type[matrix_type_size];
  matfp.read ((char*) matrix_type, matrix_type_size);

  char exp_matrix_type[] = "GPE_BSMAT_CSR";
  if (std::strcmp (matrix_type, exp_matrix_type) != 0) {
    std::cerr << "[GRAPHEE] [GPE_BSMAT_CSR] Wrong matrix format, expected GPE_BSMAT_CSR, while file is " << matrix_type << std::endl;
    exit (-1);
  }

  /* Read fileformat {BIN, SNAPPY} */
  int fileformat;
  matfp.read ((char*) &fileformat, sizeof(int));

  /* Readingg offsets of the matrix if needed */
  matfp.read ((char*) &offl, sizeof(int64_t));
  matfp.read ((char*) &offc, sizeof(int64_t));

  /* Matrix dimension */
  matfp.read ((char*) &m, sizeof(idx_t));
  matfp.read ((char*) &nnz, sizeof(idx_t));

  if ((nnz+m+1)*sizeof(idx_t) < prop.ram_limit) {
    ia = new idx_t [m+1];
    ia[0] = 0;
    ja = new idx_t [nnz];
    is_alloc = true;
  } else {
    std::cerr << "[GRAPHEE] [GPE_BSMAT_CSR] Requested size is beyond RAMLIMIT" << std::endl;
    matfp.close();
    exit (-1);
  }

  if (fileformat == BIN) {
    matfp.read ((char*) ia, (m+1)*sizeof(idx_t));
    matfp.read ((char*) ja, (nnz)*sizeof(idx_t));

  } else if (fileformat == SNAPPY) {
    bool uncomp_succeed;

    size_t ia_snappy_size;
    matfp.read ((char*) &ia_snappy_size, sizeof(size_t));

    char* ia_snappy = new char [ia_snappy_size];
    matfp.read ((char*) ia_snappy, ia_snappy_size);

    uncomp_succeed = uncompress_snappy (ia_snappy, ia_snappy_size, &ia, m);
    delete[] ia_snappy;

    if (!uncomp_succeed) {
      std::cerr << "[GRAPHEE] [GPE_BSMAT_CSR] SNAPPY uncompression of IA failed" << std::endl;

      matfp.close();
      exit (-1);
    }

    size_t ja_snappy_size;
    matfp.read ((char*) &ja_snappy_size, sizeof(size_t));

    char* ja_snappy = new char [ja_snappy_size];
    matfp.read ((char*) ja_snappy, ja_snappy_size);

    uncomp_succeed = uncompress_snappy (ja_snappy, ja_snappy_size, (char**) &ja, nnz);
    delete[] ja_snappy;

    if (!uncomp_succeed) {
      std::cerr << "[GRAPHEE] [GPE_BSMAT_CSR] SNAPPY uncompression of JA failed" << std::endl;

      matfp.close();
      exit (-1);
    }
  }
  
  matfp.close();
}

template <class idx_t>
size_t gpe_bsmat_csr<idx_t>::size () {
  return (nnz + m+1)*sizeof(idx_t);
}

template <class idx_t>
bool gpe_bsmat_csr<idx_t>::verify () {
  return (nnz == ia[m+1]);
}

} // namespace graphee

#endif // GPE_BSMAT_CSR_H
