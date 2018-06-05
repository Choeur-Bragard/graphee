#ifndef GRAPHEE_SPARSE_MATRIX_CSR_H__
#define GRAPHEE_SPARSE_MATRIX_CSR_H__

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cstring>
#include <vector>

#include "snappy/snappy.h"

#include "utils.h"
#include "properties.h"
#include "vector.h"

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
  sparseMatrixCSR () : m(0), n(0), nnz(0) {}

  sparseMatrixCSR (properties& properties, uint64_t nlines, uint64_t ncols,
      uint64_t nonzero_elems, valueT init_val = 0.) :
    props(properties), m(nlines), n(ncols), nnz(nonzero_elems) {
      bool_matrix = typeid(valueT) == typeid(bool);
      if ((nnz+m+1)*sizeof(uint64_t) < props.ram_limit) {
        if (!bool_matrix) a.resize (nnz, init_val);
        ia.resize (m+1, 0);
        ja.resize (nnz, 0);
      } else {
        print_error("Requested size is beyond \'ram_limit\'");
        exit (-1);
      }
    }

  ~sparseMatrixCSR () {}

  void fill   (uint64_t i, uint64_t j, valueT val);
  void insert (uint64_t i, uint64_t j, valueT val);
  void remove (uint64_t i, uint64_t j, valueT val);

  void save (std::string filename, int file_format = utils::BIN);
  void load (std::string filename);

  size_t size ();
  bool verify ();

  bool empty ();

  void clear ();

  template <typename vecValueT>
  vector<vecValueT>& operator* (vector<vecValueT>& rvec);
  sparseMatrixCSR<valueT>& operator* (valueT rval);

  const std::string matrixType {"sparseMatrixCSR"};

  using valueType = valueT;

  template <typename vecValueT>
  friend class vector;

private:
  std::vector<valueT> a;
  std::vector<uint64_t> ia;
  std::vector<uint64_t> ja;

  properties& props;

  bool bool_matrix;

  uint64_t m;
  uint64_t n;
  uint64_t nnz;
  uint64_t fill_id;
}; // class sparseMatrixCSR

/*! Filling the sparse matrix with sorted entries by
 * ascending lines id
 */
template <typename valueT>
void sparseMatrixCSR<valueT>::fill (uint64_t i, uint64_t j, valueT val) {
  for (uint64_t l = fill_id+1; l <= i; l++) {
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

template <typename valueT>
void sparseMatrixCSR<valueT>::save (std::string name, int fileformat) {
  std::ofstream matfp (name, std::ios_base::binary);

  size_t matrixType_size = matrixType.size();

  /* Save explicitly matrix properties */
  matfp.write (reinterpret_cast<const char*>(&matrixType_size), sizeof(size_t));
  matfp.write (reinterpret_cast<const char*>(matrixType.c_str()), matrixType_size);

  /* Save fileformat {BIN, SNAPPY} */
  matfp.write (reinterpret_cast<const char*>(&fileformat), sizeof(int));

  /* Matrix dimension */
  matfp.write (reinterpret_cast<const char*>(&m), sizeof(uint64_t));
  matfp.write (reinterpret_cast<const char*>(&nnz), sizeof(uint64_t));

  if (fileformat == utils::BIN) {
    matfp.write (reinterpret_cast<const char*>(ia.data()), ia.size()*sizeof(uint64_t));
    matfp.write (reinterpret_cast<const char*>(ja.data()), ja.size()*sizeof(uint64_t));

  } else if (fileformat == utils::SNAPPY) {
    size_t ia_snappy_size = snappy::MaxCompressedLength64 (ia.size()*sizeof(uint64_t));
    char* ia_snappy = new char [ia_snappy_size];
    snappy::RawCompress64 (reinterpret_cast<char*>(ia.data()), ia.size()*sizeof(uint64_t), ia_snappy, &ia_snappy_size);

    matfp.write (reinterpret_cast<const char*>(&ia_snappy_size), sizeof(size_t));
    matfp.write (reinterpret_cast<const char*>(ia_snappy), ia_snappy_size);
    delete[] ia_snappy;

    size_t ja_snappy_size = snappy::MaxCompressedLength64 (ja.size()*sizeof(uint64_t));
    char* ja_snappy = new char [ja_snappy_size];
    snappy::RawCompress64 (reinterpret_cast<char*>(ja.data()), ja.size()*sizeof(uint64_t), ja_snappy, &ja_snappy_size);

    matfp.write (reinterpret_cast<const char*>(&ja_snappy_size), sizeof(size_t));
    matfp.write (reinterpret_cast<const char*>(ja_snappy), ja_snappy_size);
    delete[] ja_snappy;
  }

  matfp.close();
}

template <typename valueT>
void sparseMatrixCSR<valueT>::load (std::string name) {
  this->clear ();
  std::ifstream matfp (name, std::ios_base::binary);

  /* Save explicitly matrix properties */
  size_t matrixType_size;
  matfp.read (reinterpret_cast<char*>(&matrixType_size), sizeof(size_t));

  char read_matrixType[matrixType_size+1];
  matfp.read (reinterpret_cast<char*>(read_matrixType), matrixType_size);
  read_matrixType[matrixType_size] = '\0';

  if (std::strcmp (read_matrixType, matrixType.c_str()) != 0) {
    std::ostringstream oss;
    oss << "Wrong matrix format, found \'" << read_matrixType << "\' while expecting \'"
      << matrixType << "\'";
    print_error (oss.str());
    exit (-1);
  }

  /* Read fileformat {BIN, SNAPPY} */
  int fileformat;
  matfp.read (reinterpret_cast<char*>(&fileformat), sizeof(int));

  /* Matrix dimension */
  matfp.read (reinterpret_cast<char*>(&m), sizeof(uint64_t));
  matfp.read (reinterpret_cast<char*>(&nnz), sizeof(uint64_t));

  if ((m+1)*sizeof(uint64_t) + nnz*sizeof(uint64_t) < props.ram_limit) {
    ia.resize (m+1, 0);
    ja.resize (nnz, 0);
  } else {
    print_error ("Requested size is beyond \'ram_limit\'");
    matfp.close();
    exit (-1);
  }

  if (fileformat == utils::BIN) {
    matfp.read (reinterpret_cast<char*>(ia.data()), ia.size()*sizeof(uint64_t));
    matfp.read (reinterpret_cast<char*>(ja.data()), ja.size()*sizeof(uint64_t));

  } else if (fileformat == utils::SNAPPY) {
    bool uncomp_succeed;

    size_t ia_snappy_size;
    matfp.read (reinterpret_cast<char*>(&ia_snappy_size), sizeof(size_t));

    char* ia_snappy = new char [ia_snappy_size];
    matfp.read (reinterpret_cast<char*>(ia_snappy), ia_snappy_size);

    uncomp_succeed = snappy::RawUncompress64 (ia_snappy, ia_snappy_size, reinterpret_cast<char*>(ia.data()));
    delete[] ia_snappy;

    if (!uncomp_succeed) {
      print_error("SNAPPY uncompression of IA failed");

      matfp.close();
      exit (-1);
    }

    size_t ja_snappy_size;
    matfp.read (reinterpret_cast<char*>(&ja_snappy_size), sizeof(size_t));

    char* ja_snappy = new char [ja_snappy_size];
    matfp.read (reinterpret_cast<char*>(ja_snappy), ja_snappy_size);

    uncomp_succeed = snappy::RawUncompress64 (ja_snappy, ja_snappy_size, reinterpret_cast<char*>(ja.data()));
    delete[] ja_snappy;

    if (!uncomp_succeed) {
      print_error("SNAPPY uncompression of JA failed");

      matfp.close();
      exit (-1);
    }
  }

  fill_id = m-1;

  matfp.close();
}

template <typename valueT>
size_t sparseMatrixCSR<valueT>::size () {
  return (nnz+m+1)*sizeof(uint64_t);
}

template <typename valueT>
bool sparseMatrixCSR<valueT>::verify () {
  if (fill_id < m-1) {
    for (uint64_t l = fill_id+1; l < m; l++) {
      ia[l+1] = ia[l];
    }
    fill_id = m-1;
  }

  if (nnz == ia[m]) {
    return true;
  } else {
    std::ostringstream oss;
    oss << "NNZ = " << nnz << " IA[M] = " << ia[m];
    print_warning (oss.str());
    return false;
  }
}

template <typename valueT>
void sparseMatrixCSR<valueT>::clear() {
  if (!bool_matrix) a.clear();
  ia.clear();
  ja.clear();

  m = 0;
  n = 0;
  nnz = 0;
}

template <typename valueT>
template <typename vecValueT>
vector<vecValueT>& sparseMatrixCSR<valueT>::operator* (vector<vecValueT>& rvec) {
  if (n != rvec.m) {
    std::ostringstream oss;
    oss << "Error SpMat[" << m << "x" << n << "] with Vec[" << rvec.m << "]";
    print_error (oss.str());
    exit (-1);
  }

  vector<vecValueT> res {rvec.props, m, 0.};

  if (!bool_matrix) {
    for (uint64_t i = 0; i < m; i++) {
      for (uint64_t ja_idx = ia[i]; ja_idx < ia[i+1]; ja_idx++) {
        res[i] += a[ja[ja_idx]]*rvec[ja[ja_idx]];
      }
    }
  } else {
    for (uint64_t i = 0; i < m; i++) {
      for (uint64_t ja_idx = ia[i]; ja_idx < ia[i+1]; ja_idx++) {
        res[i] += rvec[ja[ja_idx]];
      }
    }
  }
}

template <typename valueT>
sparseMatrixCSR<valueT>& sparseMatrixCSR<valueT>::operator* (valueT rval) {
  if (bool_matrix) {
    print_warning ("Scalar-matrix product applied to a \'boolean\' matrix. Nothing will be changed...");
    return *(this);
  } else {
    sparseMatrixCSR<valueT> res {props, m, n, nnz};
    for (auto& val : res) val *= rval;
    return res;
  }
}

} // namespace graphee

#endif // GRAPHEE_SPARSE_MATRIX_CSR_H__
