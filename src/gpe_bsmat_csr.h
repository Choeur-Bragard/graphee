#ifndef GPE_BSMAT_CSR_H
#define GPE_BSMAT_CSR_H

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cstring>
#include <vector>

#include <snappy.h>

#include "gpe_props.h"
#include "gpe_utils.h"

namespace graphee {

/*! \brief Boolean sparse matrix in CSR
 *         Commonly named Compressed Sparse Row
 *
 * It defines a sparse matrix filled only with elements,
 * of values \{0,1\}. This is the case of a non-weighted-edge
 * graphs.
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

  void save (std::string name, int fileformat = BIN, uint64_t offl = 0, uint64_t offc = 0);
  void load (std::string name);

  size_t size ();
  bool verify ();
  idx_t last_id {0};

private:
  gpe_props prop;

  bool is_alloc {false};

  idx_t m;
  uint64_t nnz;

  uint64_t offl;
  uint64_t offc;

  std::vector<uint64_t> ia;
  std::vector<idx_t> ja;
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

  if ((m+1)*sizeof(uint64_t)+nnz*sizeof(idx_t) < prop.ram_limit) {
    ia.resize(m+1, 0);
    ja.resize(nnz, 0);
    is_alloc = true;
  } else {
    gpe_error("Requested size is beyond \'ram_limit\'");
    exit (-1);
  }
}

/*! General destructor of the class */
template <class idx_t>
gpe_bsmat_csr<idx_t>::~gpe_bsmat_csr () {
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
void gpe_bsmat_csr<idx_t>::save (std::string name, int fileformat, uint64_t offl, uint64_t) {
  std::ofstream matfp (name, std::ios_base::binary);

  std::string matrix_type ("GPE_BSMAT_CSR");
  size_t matrix_type_size = matrix_type.size();

  /* Save explicitly matrix properties */
  matfp.write ((const char*) &matrix_type_size, sizeof(size_t));
  matfp.write ((const char*) matrix_type.c_str(), matrix_type_size);

  /* Save fileformat {BIN, SNAPPY} */
  matfp.write ((const char*) &fileformat, sizeof(int));

  /* Printing offsets of the matrix if needed */
  matfp.write ((const char*) &offl, sizeof(uint64_t));
  matfp.write ((const char*) &offc, sizeof(uint64_t));

  /* Matrix dimension */
  matfp.write ((const char*) &m, sizeof(idx_t));
  matfp.write ((const char*) &nnz, sizeof(uint64_t));

  if (fileformat == BIN) {
    matfp.write ((const char*) ia.data(), ia.size()*sizeof(uint64_t));
    matfp.write ((const char*) ja.data(), ja.size()*sizeof(idx_t));

  } else if (fileformat == SNAPPY) {
    char* ia_snappy = new char [snappy::MaxCompressedLength(ia.size()*sizeof(uint64_t))];
    size_t ia_snappy_size;
    compress_snappy ((char*)ia.data(), ia.size()*sizeof(uint64_t), ia_snappy, ia_snappy_size);

    matfp.write ((const char*) &ia_snappy_size, sizeof(size_t));
    matfp.write ((const char*) ia_snappy, ia_snappy_size);
    delete[] ia_snappy;

    char* ja_snappy = new char [snappy::MaxCompressedLength(ja.size()*sizeof(idx_t))];
    size_t ja_snappy_size;
    compress_snappy ((char*)ja.data(), (nnz)*sizeof(idx_t), ja_snappy, ja_snappy_size);

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

  char* matrix_type = new char [matrix_type_size+1];
  matfp.read ((char*) matrix_type, matrix_type_size);
  matrix_type[matrix_type_size] = '\0';

  char exp_matrix_type[] = "GPE_BSMAT_CSR";
  if (std::strcmp (matrix_type, exp_matrix_type) != 0) {
    std::ostringstream oss;
    oss << "Wrong matrix format, expected " << exp_matrix_type << " " << std::strlen(exp_matrix_type)
      << ", instead " << matrix_type << " " << std::strlen(matrix_type);
    gpe_error (oss.str());
    exit (-1);
  }

  /* Read fileformat {BIN, SNAPPY} */
  int fileformat;
  matfp.read ((char*) &fileformat, sizeof(int));

  /* Readingg offsets of the matrix if needed */
  matfp.read ((char*) &offl, sizeof(uint64_t));
  matfp.read ((char*) &offc, sizeof(uint64_t));

  /* Matrix dimension */
  matfp.read ((char*) &m, sizeof(idx_t));
  matfp.read ((char*) &nnz, sizeof(uint64_t));

  if ((m+1)*sizeof(uint64_t) + nnz*sizeof(idx_t) < prop.ram_limit) {
    ia.resize (m+1, 0);
    ja.resize (nnz, 0);
    is_alloc = true;
  } else {
    gpe_error ("Requested size is beyond \'RAMLIMIT\'");
    matfp.close();
    exit (-1);
  }

  if (fileformat == BIN) {
    matfp.read ((char*) ia.data(), ia.size()*sizeof(uint64_t));
    matfp.read ((char*) ja.data(), ja.size()*sizeof(idx_t));

  } else if (fileformat == SNAPPY) {
    bool uncomp_succeed;

    size_t ia_snappy_size;
    matfp.read ((char*) &ia_snappy_size, sizeof(size_t));

    char* ia_snappy = new char [ia_snappy_size];
    matfp.read ((char*) ia_snappy, ia_snappy_size);

    uncomp_succeed = uncompress_snappy (ia_snappy, ia_snappy_size, (char*)ia.data(), m+1);
    delete[] ia_snappy;

    if (!uncomp_succeed) {
      gpe_error("SNAPPY uncompression of IA failed");

      matfp.close();
      exit (-1);
    }

    size_t ja_snappy_size;
    matfp.read ((char*) &ja_snappy_size, sizeof(size_t));

    char* ja_snappy = new char [ja_snappy_size];
    matfp.read ((char*) ja_snappy, ja_snappy_size);

    uncomp_succeed = uncompress_snappy (ja_snappy, ja_snappy_size, (char*)ja.data(), nnz);
    delete[] ja_snappy;

    if (!uncomp_succeed) {
      gpe_error("SNAPPY uncompression of JA failed");

      matfp.close();
      exit (-1);
    }
  }

  last_id = m;

  matfp.close();
}

template <class idx_t>
size_t gpe_bsmat_csr<idx_t>::size () {
  return (m+1)*sizeof(uint64_t) + nnz*sizeof(idx_t);
}

template <class idx_t>
bool gpe_bsmat_csr<idx_t>::verify () {
  if (last_id < m) {
    for (idx_t l = last_id+1; l <= m; l++) {
      ia[l+1] = ia[l];
    }
    last_id = m;
  }

  if (nnz == ia[m]) {
    return true;
  } else {
    std::ostringstream oss;
    oss << "Wrong matrix : M = " << m << "; NNZ = " << nnz << "; IA[M] = " << ia[m];
    gpe_error (oss.str());
    return false;
  }
}

} // namespace graphee

#endif // GPE_BSMAT_CSR_H
