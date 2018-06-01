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
  sparseMatrixCSR (properties& properties, uint64_t nlines, uint64_t ncols,
      uint64_t nonzero_elems, valueT init_val = 0.);
  ~sparseMatrixCSR ();

  void fill   (uint64_t i, uint64_t j, valueT val);
  void insert (uint64_t i, uint64_t j, valueT val);
  void remove (uint64_t i, uint64_t j, valueT val);

  size_t size () const;
  bool verify () const;

  void clear ();

  std::string& get_matrix_properties ();

  template <typename vecValueT>
  vector<vecValueT>& operator* (sparseMatrixCSR<valueT>& lmat, vector<vecValueT>& rvec);

  sparseMatrixCSR<valueT>& operator* (valueT lval, sparseMatrixCSR<valueT>& rmat);
  sparseMatrixCSR<valueT>& operator* (sparseMatrixCSR<valueT>& lmat, valueT rval);

  const std::string matrixType {"sparseMatrixCSR"};

  friend vector<valueT>;

private:
  std::vector<valueT> a;
  std::vector<uint64_t> ia;
  std::vector<uint64_t> ja;

  properties&& props;

  bool bool_matrix;

  uint64_t m;
  uint64_t n;
  uint64_t nnz;
  uint64_t fill_id;
}; // class sparseMatrixCSR

/*! Empty constructor */
template <typename valueT>
sparseMatrixCSR<valueT>::sparseMatrixCSR () {
  bool_matrix {std::typeid(valueT) != std::typeid(bool)};
  m {0};
  n {0};
  nnz {0};
  fill_id {0};
}

/*! General constructor of class */
template <typename valueT>
sparseMatrixCSR<valueT>::sparseMatrixCSR (properties& properties, uint64_t nlines, uint64_t nonzero_elems, valueT init_val);
  props {properties};
  m {nlines};
  n {ncols};
  nnz {nonzero_elems};

  if ((nzz+m+1)*sizeof(uint64_t) < props.ram_limit) {
    if (!bool_matrix) a.resize (nnz, init_val);
    ia.resize (m+1, 0);
    ja.resize (nnz, 0);
  } else {
    print_error("Requested size is beyond \'ram_limit\'");
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


template <typename valueT>
size_t sparseMatrixCSR<valueT>::size () const {
  return (nnz+m+1)*sizeof(uint64_t);
}

template <typename valueT>
bool sparseMatrixCSR<valueT>::verify () const {
  if (fill_id < m) {
    for (uint64_t l = fill_id+1; l <= m; l++) {
      ia[l+1] = ia[l];
    }
    fill_id = m;
  }

  if (nnz == ia[m]) {
    return true;
  } else {
    std::ostringstream oss;
    oss << "NNZ = " << nnz << " IA[M+1] = " << ia[m+1];
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
vector<valueT>& operator* (sparseMatrixCSR<valueT>& lmat, vector<valueT>& rvec) {
  if (lmat.n != rvec.m) {
    std::ostringstream oss;
    oss << "Error SpMat[" << lmat.m << "x" << lmat.n << "] with Vec[" << rvec.m << "]";
    print_error (oss.str());
    exit (-1);
  }

  vector<valueT> res {rvec.props, lmat.m, 0.};

  if (!bool_matrix) {
    for (uint64_t i = 0; i < lmat.m; i++) {
      for (uint64_t ja_idx = ia[i]; ja_idx < ia[i+1]; ja_idx++) {
        res[i] += lmat.a[ lmat.ja[ja_idx] ]*rvec[j];
      }
    }
  } else {
    for (uint64_t i = 0; i < lmat.m; i++) {
      for (uint64_t ja_idx = ia[i]; ja_idx < ia[i+1]; ja_idx++) {
        res[i] += rvec[j];
      }
    }
  }
}

template <typename valueT>
sparseMatrixCSR<valueT>& operator* (valueT lval, sparseMatrixCSR<valueT>& rmat) {
  if (bool_matrix) {
    print_warning ("Scalar-matrix product applied to a \'boolean\' matrix. Nothing will be changed...");
    return rmat;
  } else {
    sparseMatrixCSR<valueT> res {rmat};
    for (auto& val : res) val *= lval;
    return res;
  }
}

template <typename valueT>
sparseMatrixCSR<valueT>& operator* (sparseMatrixCSR<valueT>& lmat, valueT rval) {
  if (bool_matrix) {
    print_warning ("Scalar-matrix product applied to a \'boolean\' matrix. Nothing will be changed...");
    return lmat;
  } else {
    sparseMatrixCSR<valueT> res {lmat};
    for (auto& val : res) val *= rval;
    return res;
  }
}

} // namespace graphee

#endif // sparseMatrixCSR_H
