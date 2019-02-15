#ifndef GRAPHEE_SPARSE_MATRIX_CSR_HPP__
#define GRAPHEE_SPARSE_MATRIX_CSR_HPP__

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cstring>
#include <vector>
#include <typeinfo>

#include "snappy/snappy.h"

#include "utils.hpp"
#include "properties.hpp"

#include "sparse_bmatrix_csr.hpp"

namespace graphee
{

/*! \brief Sparse matrix in CSR format
 *         Commonly named Compressed Sparse Row
 *
 * The matrix can handle either `double`, `float` and `int`
 * values.
 */

template <typename ValueT>
class SparseMatrixCSR : public SparseBMatrixCSR
{
public:
  SparseMatrixCSR(Properties *properties) : SparseBMatrixCSR(properties)
  {
    if (typeid(ValueT) == typeid(bool))
    {
      print_error ("Do not use \'SparseMatrixCSR<bool>\', instead use \'SparseBMatrixCSR\'");
      exit(-1);
    }
  }

  SparseMatrixCSR(Properties *properties, uint64_t nlines, uint64_t ncols,
                  uint64_t nonzero_elems, ValueT init_val = 0.) :
    SparseBMatrixCSR(properties, nlines, ncols, nonzero_elems)
  {
    if (typeid(ValueT) == typeid(bool))
    {
      print_error ("Do not use \'SparseMatrixCSR<bool>\', instead use \'SparseBMatrixCSR\'");
      exit(-1);
    }

    if ((nnz + m + 1) * sizeof(uint64_t) + nnz * sizeof(ValueT) < props->ram_limit)
    {
      a.resize(nnz, init_val);
    }
    else
    {
      print_error("Requested size is beyond \'ram_limit\'");
      exit(-1);
    }
  }

  SparseMatrixCSR(SparseMatrixCSR<ValueT> &&mat) : SparseBMatrixCSR(std::move(mat)), a(std::move(mat.a)) {}

  ~SparseMatrixCSR()
  {
    a.clear();
  }

  void fill(uint64_t i, uint64_t j, ValueT val);
  void insert(uint64_t i, uint64_t j, ValueT val);
  void remove(uint64_t i, uint64_t j, ValueT val);

  void save(std::string filename, int file_format = Utils::BIN);
  void load(std::string filename);

  size_t size();

  bool empty();

  void clear();

  template <typename vecValueT>
  Vector<vecValueT> &operator*(Vector<vecValueT> &rvec);
  SparseMatrixCSR<ValueT> &operator*(ValueT rval);

  Vector<float> columns_sum();

  const std::string matrix_typename{"SparseMatrixCSR"};

  using ValueType = ValueT;

private:
  std::vector<ValueT> a;
}; // class SparseMatrixCSR

/*! Filling the sparse matrix with sorted entries by
 * ascending lines id
 */
template <typename ValueT>
void SparseMatrixCSR<ValueT>::fill(uint64_t i, uint64_t j, ValueT val)
{
  SparseBMatrixCSR::fill(i, j);
  a[ia[i + 1]] = val;
}

/*! Inserting element in a CSR matrix */
template <typename ValueT>
void SparseMatrixCSR<ValueT>::insert(uint64_t i, uint64_t j, ValueT val)
{
}

/*! Remove element of the CSR matrix*/
template <typename ValueT>
void SparseMatrixCSR<ValueT>::remove(uint64_t i, uint64_t j, ValueT val)
{
}

template <typename ValueT>
void SparseMatrixCSR<ValueT>::save(std::string name, int fileformat)
{
  std::ofstream matfp(name, std::ios_base::binary);

  size_t matrix_typename_size = matrix_typename.size();

  /* Save explicitly matrix properties */
  matfp.write(reinterpret_cast<const char *>(&matrix_typename_size), sizeof(size_t));
  matfp.write(reinterpret_cast<const char *>(matrix_typename.c_str()), matrix_typename_size);

  /* Save fileformat {BIN, SNAPPY} */
  matfp.write(reinterpret_cast<const char *>(&fileformat), sizeof(int));

  /* Matrix dimension */
  matfp.write(reinterpret_cast<const char *>(&m), sizeof(uint64_t));
  matfp.write(reinterpret_cast<const char *>(&nnz), sizeof(uint64_t));

  if (fileformat == Utils::BIN)
  {
    matfp.write(reinterpret_cast<const char *>(a.data()), a.size() * sizeof(ValueT));
    matfp.write(reinterpret_cast<const char *>(ia.data()), ia.size() * sizeof(uint64_t));
    matfp.write(reinterpret_cast<const char *>(ja.data()), ja.size() * sizeof(uint64_t));
  }
  else if (fileformat == Utils::SNAPPY)
  {
    size_t a_snappy_size = snappy::MaxCompressedLength64(a.size() * sizeof(ValueT));
    char *a_snappy = new char[a_snappy_size];
    snappy::RawCompress64(reinterpret_cast<char *>(a.data()), a.size() * sizeof(ValueT), a_snappy, &a_snappy_size);

    matfp.write(reinterpret_cast<const char *>(&a_snappy_size), sizeof(size_t));
    matfp.write(reinterpret_cast<const char *>(a_snappy), a_snappy_size);
    delete[] a_snappy;

    size_t ia_snappy_size = snappy::MaxCompressedLength64(ia.size() * sizeof(uint64_t));
    char *ia_snappy = new char[ia_snappy_size];
    snappy::RawCompress64(reinterpret_cast<char *>(ia.data()), ia.size() * sizeof(uint64_t), ia_snappy, &ia_snappy_size);

    matfp.write(reinterpret_cast<const char *>(&ia_snappy_size), sizeof(size_t));
    matfp.write(reinterpret_cast<const char *>(ia_snappy), ia_snappy_size);
    delete[] ia_snappy;

    size_t ja_snappy_size = snappy::MaxCompressedLength64(ja.size() * sizeof(uint64_t));
    char *ja_snappy = new char[ja_snappy_size];
    snappy::RawCompress64(reinterpret_cast<char *>(ja.data()), ja.size() * sizeof(uint64_t), ja_snappy, &ja_snappy_size);

    matfp.write(reinterpret_cast<const char *>(&ja_snappy_size), sizeof(size_t));
    matfp.write(reinterpret_cast<const char *>(ja_snappy), ja_snappy_size);
    delete[] ja_snappy;
  }

  matfp.close();
}

template <typename ValueT>
void SparseMatrixCSR<ValueT>::load(std::string name)
{
  this->clear();
  std::ifstream matfp(name, std::ios_base::binary);

  /* Save explicitly matrix properties */
  size_t matrix_typename_size;
  matfp.read(reinterpret_cast<char *>(&matrix_typename_size), sizeof(size_t));

  char read_matrix_typename[matrix_typename_size + 1];
  matfp.read(reinterpret_cast<char *>(read_matrix_typename), matrix_typename_size);
  read_matrix_typename[matrix_typename_size] = '\0';

  if (std::strcmp(read_matrix_typename, matrix_typename.c_str()) != 0)
  {
    std::ostringstream oss;
    oss << "Wrong matrix format, found \'" << read_matrix_typename << "\' while expecting \'"
        << matrix_typename << "\'";
    print_error(oss.str());
    exit(-1);
  }

  /* Read fileformat {BIN, SNAPPY} */
  int fileformat;
  matfp.read(reinterpret_cast<char *>(&fileformat), sizeof(int));

  /* Matrix dimension */
  matfp.read(reinterpret_cast<char *>(&m), sizeof(uint64_t));
  matfp.read(reinterpret_cast<char *>(&nnz), sizeof(uint64_t));

  if ((nnz + m + 1) * sizeof(uint64_t) + nnz * sizeof(ValueT) < props->ram_limit)
  {
    a.resize(nnz, 0.);
    ia.resize(m + 1, 0);
    ja.resize(nnz, 0);
  }
  else
  {
    print_error("Requested size is beyond \'ram_limit\'");
    matfp.close();
    exit(-1);
  }

  if (fileformat == Utils::BIN)
  {
    matfp.read(reinterpret_cast<char *>(a.data()), a.size() * sizeof(ValueT));
    matfp.read(reinterpret_cast<char *>(ia.data()), ia.size() * sizeof(uint64_t));
    matfp.read(reinterpret_cast<char *>(ja.data()), ja.size() * sizeof(uint64_t));
  }
  else if (fileformat == Utils::SNAPPY)
  {
    bool uncomp_succeed;

    size_t a_snappy_size;
    matfp.read(reinterpret_cast<char *>(&a_snappy_size), sizeof(size_t));

    char *a_snappy = new char[a_snappy_size];
    matfp.read(reinterpret_cast<char *>(a_snappy), a_snappy_size);

    uncomp_succeed = snappy::RawUncompress64(a_snappy, a_snappy_size, reinterpret_cast<char *>(a.data()));
    delete[] a_snappy;

    if (!uncomp_succeed)
    {
      print_error("SNAPPY uncompression of A failed");

      matfp.close();
      exit(-1);
    }

    size_t ia_snappy_size;
    matfp.read(reinterpret_cast<char *>(&ia_snappy_size), sizeof(size_t));

    char *ia_snappy = new char[ia_snappy_size];
    matfp.read(reinterpret_cast<char *>(ia_snappy), ia_snappy_size);

    uncomp_succeed = snappy::RawUncompress64(ia_snappy, ia_snappy_size, reinterpret_cast<char *>(ia.data()));
    delete[] ia_snappy;

    if (!uncomp_succeed)
    {
      print_error("SNAPPY uncompression of IA failed");

      matfp.close();
      exit(-1);
    }

    size_t ja_snappy_size;
    matfp.read(reinterpret_cast<char *>(&ja_snappy_size), sizeof(size_t));

    char *ja_snappy = new char[ja_snappy_size];
    matfp.read(reinterpret_cast<char *>(ja_snappy), ja_snappy_size);

    uncomp_succeed = snappy::RawUncompress64(ja_snappy, ja_snappy_size, reinterpret_cast<char *>(ja.data()));
    delete[] ja_snappy;

    if (!uncomp_succeed)
    {
      print_error("SNAPPY uncompression of JA failed");

      matfp.close();
      exit(-1);
    }
  }

  fill_id = m - 1;

  matfp.close();
}

template <typename ValueT>
size_t SparseMatrixCSR<ValueT>::size()
{
  return (nnz + m + 1) * sizeof(uint64_t) + nnz * sizeof(ValueT);
}


template <typename ValueT>
void SparseMatrixCSR<ValueT>::clear()
{
  SparseBMatrixCSR::clear();
  a.clear();
}

template <typename ValueT>
template <typename vecValueT>
Vector<vecValueT> &SparseMatrixCSR<ValueT>::operator*(Vector<vecValueT> &rvec)
{
  if (n != rvec.get_lines())
  {
    std::ostringstream oss;
    oss << "Error SpMat[" << m << "x" << n << "] with Vec[" << rvec.get_lines() << "]";
    print_error(oss.str());
    exit(-1);
  }

  Vector<vecValueT> res(rvec.get_properties(), m, 0.);

#pragma omp parallel for num_threads(props->nthreads)
  for (uint64_t i = 0; i < m; i++)
  {
    for (uint64_t ja_idx = ia[i]; ja_idx < ia[i + 1]; ja_idx++)
    {
      res[i] += a[ja[ja_idx]] * rvec[ja[ja_idx]];
    }
  }
}

template <typename ValueT>
Vector<float> SparseMatrixCSR<ValueT>::columns_sum(){

  Vector<float> res(props, m, 0.);

#pragma omp parallel for num_threads (props->nthreads)
  for(uint64_t i = 0; i < ja.size(); i++){
    res[ja[i]]+=1;
  }

  return res;
}


template <typename ValueT>
SparseMatrixCSR<ValueT> &SparseMatrixCSR<ValueT>::operator*(ValueT rval)
{
  for (uint64_t i = 0; i < nnz; i++)
    a[i] *= rval;

  return *this;
}

} // namespace graphee

#endif // GRAPHEE_SPARSE_MATRIX_CSR_HPP__
