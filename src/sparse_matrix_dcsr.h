#ifndef GRAPHEE_SPARSE_MATRIX_DCSR_H__ 
#define GRAPHEE_SPARSE_MATRIX_DCSR_H__ 

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <cstring>
#include <limits>

#include "snappy/snappy.h"

#include "graphee.h"

namespace graphee {

/*! \brief Boolean sparse matrix in DCSR
 *         Commonly named Dynamic Compressed Sparse Row
 *
 *  It defines a sparse matrix with boolean values either 0 or 1,
 *  thus no array `a` is allocated. This could be used for 
 *  non-wighted-edge graphs.
 */

template <typename idx_t>
class gpe_bsmat_dcsr {
public:
  enum {BIN, SNAPPY};

  gpe_bsmat_dcsr (gpe_props& i_prop);
  gpe_bsmat_dcsr (gpe_props& i_prop, idx_t i_m, idx_t nnz, uint64_t i_pitch = 2, uint64_t i_alpha = 200);
  ~gpe_bsmat_dcsr ();

  void sorted_fill (idx_t i, idx_t j);

  void set_offsets (uint64_t offl, uint64_t offc);

  void insert (idx_t i, idx_t j);
  void remove (idx_t i, idx_t j);
  void defrag ();

  void save (std::string name, int fileformat = BIN, uint64_t offl = 0, uint64_t offc = 0);
  void load (std::string name);

  size_t size ();
  bool verify ();

private:
  gpe_props prop;

  bool is_alloc {false};

  idx_t last_id {0};

  std::vector<uint64_t> ia;
  std::vector<idx_t> ja;

  uint64_t data_ptr {0};
  uint64_t pitch {0};
  uint64_t alpha {0};

  idx_t m;
  uint64_t nnz;

  uint64_t offl;
  uint64_t offc;
}; // gpe_bsmat_dcsr

template <typename idx_t>
gpe_bsmat_dcsr<idx_t>::gpe_bsmat_dcsr (gpe_props& i_prop) {
  if (i_prop.window > std::numeric_limits<idx_t>::max()) {
    gpe_error ("\'idx_t\' has a lower limit than the \'window\' value");
    exit (-1);
  }
  prop = i_prop;
}

template <typename idx_t>
gpe_bsmat_dcsr<idx_t>::gpe_bsmat_dcsr (gpe_props& i_prop, idx_t i_m, idx_t i_nnz, size_t i_pitch, size_t i_alpha) {
  if (i_prop.window > std::numeric_limits<idx_t>::max()) {
    gpe_error ("\'idx_t\' has a lower limit than the \'window\' value");
    exit (-1);
  }

  prop = i_prop;
  m = i_m;
  nnz = i_nnz;
  pitch = i_pitch;
  alpha = i_alpha;

  ia.resize (2*m*pitch, 0);
  ja.resize (nnz+alpha, 0);
}

template <typename idx_t>
void gpe_bsmat_dcsr<idx_t>::sorted_fill (idx_t i, idx_t j) {
  for (idx_t l = last_id; l <= i; l++) {
    ia[2*(l+1)*pitch]   = ia[2*l*pitch+1];
    ia[2*(l+1)*pitch+1] = ia[2*l*pitch+1];
  }

  ja[ia[2*i*pitch+1]] = j;
  ia[2*i*pitch+1]++;
  data_ptr++;

  last_id = i;
}

template <typename idx_t>
void gpe_bsmat_dcsr<idx_t>::set_offsets (uint64_t i_offl, uint64_t i_offc) {
  offl = i_offl;
  offc = i_offc;
}

/*! Inserting element in a CSR matrix */
template <typename idx_t>
void gpe_bsmat_dcsr<idx_t>::insert (idx_t i, idx_t j) {
}

/*! Remove element of the CSR matrix*/
template <typename idx_t>
void gpe_bsmat_dcsr<idx_t>::remove (idx_t i, idx_t j) {
}

template <typename idx_t>
void gpe_bsmat_dcsr<idx_t>::defrag () {
}

/*! Save the matrix to a file
 *
 * One can determine the fileformat avail BIN or
 * SNAPPY for fast-light compression.
 */
template <typename idx_t>
void gpe_bsmat_dcsr<idx_t>::save (std::string name, int fileformat, uint64_t offl, uint64_t offc) {
  std::ofstream matfp (name, std::ios_base::binary);

  std::string matrix_type ("GPE_BSMAT_DCSR");
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
  matfp.write ((const char*) &nnz, sizeof(uint64_t));
  matfp.write ((const char*) &pitch, sizeof(uint64_t));
  matfp.write ((const char*) &alpha, sizeof(uint64_t));

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
    compress_snappy ((char*)ja.data(), ja.size()*sizeof(idx_t), ja_snappy, ja_snappy_size);

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
template <typename idx_t>
void gpe_bsmat_dcsr<idx_t>::load (std::string name) {
  std::ifstream matfp (name, std::ios_base::binary);

  /* Save explicitly matrix properties */
  size_t matrix_type_size;
  matfp.read ((char*) &matrix_type_size, sizeof(size_t));

  char matrix_type[matrix_type_size];
  matfp.read ((char*) matrix_type, matrix_type_size);

  char exp_matrix_type[] = "GPE_BSMAT_DCSR";
  if (std::strcmp (matrix_type, exp_matrix_type) != 0) {
    gpe_error ("Wrong matrix format, expected GPE_BSMAT_DCSR");
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
  matfp.read ((char*) &nnz, sizeof(uint64_t));
  matfp.read ((char*) &pitch, sizeof(uint64_t));
  matfp.read ((char*) &alpha, sizeof(uint64_t));

  if ((2*m*pitch)*sizeof(uint64_t)+(nnz+alpha)*sizeof(idx_t) < prop.ram_limit) {
    ia.resize (2*m*pitch, 0);
    ja.resize (nnz+alpha, 0);
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

    uncomp_succeed = uncompress_snappy (ia_snappy, ia_snappy_size, (char*)ia.data(), 2*m*pitch);
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

    uncomp_succeed = uncompress_snappy (ja_snappy, ja_snappy_size, (char*)ja.data(), nnz+alpha);
    delete[] ja_snappy;

    if (!uncomp_succeed) {
      gpe_error("SNAPPY uncompression of JA failed");

      matfp.close();
      exit (-1);
    }
  }
  
  matfp.close();
}

template <typename idx_t>
size_t gpe_bsmat_dcsr<idx_t>::size () {
  return (2*m*pitch)*sizeof(uint64_t)+(nnz+alpha)*sizeof(idx_t);
}

template <typename idx_t>
bool gpe_bsmat_dcsr<idx_t>::verify () {
  return (nnz == ia[2*m*pitch]);
}
} // graphee

#endif // GRAPHEE_SPARSE_MATRIX_DCSR_H__
