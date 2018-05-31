#ifndef SPARSE_MATRIX_CSR_H
#define SPARSE_MATRIX_CSR_H

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cstring>
#include <vector>

#include <snappy.h>

#include "graphee.h"

namespace graphee {

/*! \brief Sparse matrix in CSR format
 *         Commonly named Compressed Sparse Row
 *
 * The matrix can handle either `double`, `float` and `int` 
 * values.
 * A special case is when these matrices save 
 * `bool` values. In such cases no vector `a` is allocated.
 */

template <typename valueT>
class sparseMatrixCSR {
public:
  typedef valueT valueType;

  sparseMatrixCSR ();
  sparseMatrixCSR (properties& properties, uint64_t nlines, uint64_t nonzero_elems);
  ~sparseMatrixCSR ();

  void fill   (uint64_t i, uint64_t j, valueT val);
  void insert (uint64_t i, uint64_t j, valueT val);
  void remove (uint64_t i, uint64_t j, valueT val);

  size_t size () const;
  bool verify () const;

  void clear ();

  const std::string matrixType {"sparseMatrixCSR"};

private:
  std::vector<valueT> a;
  std::vector<uint64_t> ia;
  std::vector<uint64_t> ja;

  properties&& props;

  bool is_bool_sparse_matrix;

  uint64_t m;
  uint64_t nnz;
  uint64_t fill_id;
}; // class sparseMatrixCSR

/*! Empty constructor */
template <typename valueT>
sparseMatrixCSR<valueT>::sparseMatrixCSR () {
  bool_matrix {std::typeid(valueT) != std::typeid(bool)};
  m {0};
  nnz {0};
  fill_id {0};
}

/*! General constructor of class */
template <typename valueT>
sparseMatrixCSR<valueT>::sparseMatrixCSR (properties& properties, uint64_t nlines, uint64_t nonzero_elems);
  props {properties};
  m {nlines};
  nnz {nonzero_elems};

  if ((nzz+m+1)*sizeof(uint64_t) < props.ram_limit) {
    if (!bool_matrix) a.resize (nnz, 0.);
    ia.resize (m+1, 0);
    ja.resize (nnz, 0);
  } else {
    gpe_error("Requested size is beyond \'ram_limit\'");
    exit (-1);
  }
}

/*! General destructor of the class */
template <typename valueT>
sparseMatrixCSR<valueT>::~sparseMatrixCSR () {
  delete a;
  delete ia;
  delete ja;
}

/*! Filling the sparse matrix with sorted entries by
 * ascending lines id
 */
template <typename valueT>
void sparseMatrixCSR<valueT>::fill (uint64_t i, uint64_t j, valueT val) {
  for (uint64_t l = last_id+1; l <= i; l++) {
    ia[l+1] = ia[l];
  }

  if (!bool_matrix) a[ia[i+1]] = val;
  ja[ia[i+1]] = j;

  ia[i+1]++;

  fill_id = i;
}

/*! Inserting element in a CSR matrix */
template <typename valueT>
void sparseMatrixCSR<valueT>::insert (uint64_t i, uint64_t j, valueT val) {
}

/*! Remove element of the CSR matrix*/
template <typename valueT>
void sparseMatrixCSR<valueT>::remove (uint64_t i, uint64_t j, valueT val) {
}

/*! To be pushed within 'diskSparseMatrix' */
template <typename valueT>
void sparseMatrixCSR<valueT>::save (std::string name, int fileformat, uint64_t offl, uint64_t) {
  std::ofstream matfp (name, std::ios_base::binary);

  size_t matrix_type_size = matrix_type.size();

  /* Save explicitly matrix properties */
  matfp.write (reinterpret_cast<const char*>(&matrix_type_size), sizeof(size_t));
  matfp.write (reinterpret_cast<const char*>(matrix_type.c_str()), matrix_type_size);

  /* Save fileformat {BIN, SNAPPY} */
  matfp.write (reinterpret_cast<const char*>(&fileformat), sizeof(int));

  /* Printing offsets of the matrix if needed */
  matfp.write (reinterpret_cast<const char*>(&offl), sizeof(uint64_t));
  matfp.write (reinterpret_cast<const char*>(&offc), sizeof(uint64_t));

  /* Matrix dimension */
  matfp.write (reinterpret_cast<const char*>(&m), sizeof(uint64_t));
  matfp.write (reinterpret_cast<const char*>(&nnz), sizeof(uint64_t));

  if (fileformat == BIN) {
    matfp.write (reinterpret_cast<const char*>(ia.data()), ia.size()*sizeof(uint64_t));
    matfp.write (reinterpret_cast<const char*>(ja.data()), ja.size()*sizeof(uint64_t));

  } else if (fileformat == SNAPPY) {
    size_t ia_snappy_size = max_compress_size (ia.size()*sizeof(uint64_t));
    char* ia_snappy = new char [ia_snappy_size];
    compress_snappy (reinterpret_cast<char*>(ia.data()), ia.size()*sizeof(uint64_t), ia_snappy, ia_snappy_size);

    matfp.write (reinterpret_cast<const char*>(&ia_snappy_size), sizeof(size_t));
    matfp.write (reinterpret_cast<const char*>(ia_snappy), ia_snappy_size);
    delete[] ia_snappy;

    size_t ja_snappy_size = max_compress_size (ja.size()*sizeof(uint64_t));
    char* ja_snappy = new char [ja_snappy_size];
    compress_snappy (reinterpret_cast<char*>(ja.data()), (nnz)*sizeof(uint64_t), ja_snappy, ja_snappy_size);

    matfp.write (reinterpret_cast<const char*>(&ja_snappy_size), sizeof(size_t));
    matfp.write (reinterpret_cast<const char*>(ja_snappy), ja_snappy_size);
    delete[] ja_snappy;
  }

  matfp.close();
}

/*! To be pushed within 'diskSparseMatrix' */
template <typename valueT>
void sparseMatrixCSR<valueT>::load (std::string name) {
  this->clear ();
  std::ifstream matfp (name, std::ios_base::binary);

  /* Save explicitly matrix properties */
  size_t matrix_type_size;
  matfp.read (reinterpret_cast<char*>(&matrix_type_size), sizeof(size_t));

  char read_matrix_type[matrix_type_size+1];
  matfp.read (reinterpret_cast<char*>(read_matrix_type), matrix_type_size);
  read_matrix_type[matrix_type_size] = '\0';

  if (std::strcmp (read_matrix_type, matrix_type.c_str()) != 0) {
    err.str("");
    err << "Wrong matrix format, found \'" << read_matrix_type << "\' while expecting \'"
      << matrix_type << "\'";
    gpe_error (err.str());
    exit (-1);
  }

  /* Read fileformat {BIN, SNAPPY} */
  int fileformat;
  matfp.read (reinterpret_cast<char*>(&fileformat), sizeof(int));

  /* Readingg offsets of the matrix if needed */
  matfp.read (reinterpret_cast<char*>(&offl), sizeof(uint64_t));
  matfp.read (reinterpret_cast<char*>(&offc), sizeof(uint64_t));

  /* Matrix dimension */
  matfp.read (reinterpret_cast<char*>(&m), sizeof(uint64_t));
  matfp.read (reinterpret_cast<char*>(&nnz), sizeof(uint64_t));

  if ((m+1)*sizeof(uint64_t) + nnz*sizeof(uint64_t) < prop.ram_limit) {
    ia.resize (m+1, 0);
    ja.resize (nnz, 0);
    is_alloc = true;
  } else {
    gpe_error ("Requested size is beyond \'ram_limit\'");
    matfp.close();
    exit (-1);
  }

  if (fileformat == BIN) {
    matfp.read (reinterpret_cast<char*>(ia.data()), ia.size()*sizeof(uint64_t));
    matfp.read (reinterpret_cast<char*>(ja.data()), ja.size()*sizeof(uint64_t));

  } else if (fileformat == SNAPPY) {
    bool uncomp_succeed;

    size_t ia_snappy_size;
    matfp.read (reinterpret_cast<char*>(&ia_snappy_size), sizeof(size_t));

    char* ia_snappy = new char [ia_snappy_size];
    matfp.read (reinterpret_cast<char*>(ia_snappy), ia_snappy_size);

    uncomp_succeed = uncompress_snappy (ia_snappy, ia_snappy_size, reinterpret_cast<char*>(ia.data()), (m+1)*sizeof(uint64_t));
    delete[] ia_snappy;

    if (!uncomp_succeed) {
      gpe_error("SNAPPY uncompression of IA failed");

      matfp.close();
      exit (-1);
    }

    size_t ja_snappy_size;
    matfp.read (reinterpret_cast<char*>(&ja_snappy_size), sizeof(size_t));

    char* ja_snappy = new char [ja_snappy_size];
    matfp.read (reinterpret_cast<char*>(ja_snappy), ja_snappy_size);

    uncomp_succeed = uncompress_snappy (ja_snappy, ja_snappy_size, reinterpret_cast<char*>(ja.data()), nnz*sizeof(uint64_t));
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

template <typename valueT>
size_t sparseMatrixCSR<valueT>::size () {
  return (nnz+m+1)*sizeof(uint64_t);
}

template <typename valueT>
bool sparseMatrixCSR<valueT>::verify () {
  if (fill_id < m) {
    for (uint64_t l = fill_id+1; l <= m; l++) {
      ia[l+1] = ia[l];
    }
    fill_id = m;
  }

  if (nnz == ia[m]) {
    return true;
  } else {
    wrn.str("");
    wrn << "NNZ = " << nnz << " IA[M+1] = " << ia[m+1];
    gpe_warning (wrn.str());
    return false;
  }
}

template <typename valueT>
void sparseMatrixCSR<valueT>::clear() {
  if (!bool_matrix) a.clear();
  ia.clear();
  ja.clear();

  m = 0;
  nnz = 0;
}

} // namespace graphee

#endif // sparseMatrixCSR_H
