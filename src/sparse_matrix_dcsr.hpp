#ifndef GRAPHEE_SPARSE_MATRIX_DCSR_HPP__
#define GRAPHEE_SPARSE_MATRIX_DCSR_HPP__

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

#include "sparse_bmatrix_dcsr.hpp"

namespace graphee
{

/*! \brief Sparse matrix in Dynamic CSR format
*         Commonly named Dynamic Compressed Sparse Row:
*         https://thomas.gilray.org/pdf/dynamic-csr.pdf
*
* The matrix can handle either `double`, `float` and `int`
* values.
*/

template <typename ValueT>
class SparseMatrixDCSR : public SparseBMatrixDCSR
{
public:
  SparseMatrixDCSR(Properties *properties) : SparseBMatrixDCSR(properties)
  {
    if (typeid(ValueT) == typeid(bool))
    {
      print_error ("Do not use \'SparseMatrixDCSR<bool>\', instead use \'SparseBMatrixDCSR\'");
      exit(-1);
    }
  }

  SparseMatrixCSR(Properties *properties, uint64_t nlines, uint64_t ncols,
                  uint64_t nonzero_elems, uint64_t pitch, uint64_t alpha,
                  ValueT init_val = 0.) :
    SparseBMatrixCSR(properties, nlines, ncols, nonzero_elems, pitch, alpha)
  {
    if (typeid(ValueT) == typeid(bool))
    {
      print_error ("Do not use \'SparseMatrixDCSR<bool>\', instead use \'SparseBMatrixDCSR\'");
      exit(-1);
    }

    if ((nnz + alpha + 2 * m * pitch) * sizeof(uint64_t)
        + (nnz + alpha) * sizeof(ValueT) < props->ram_limit)
    {
      a.resize(nnz + alpha, init_val);
    }
    else
    {
      print_error("Requested size is beyond \'ram_limit\'");
      exit(-1);
    }
  }

  SparseMatrixDCSR(SparseMatrixDCSR<ValueT> &&mat) : SparseBMatrixDCSR(std::move(mat)), a(std::move(mat.a)) {}

  ~SparseMatrixCSR()
  {
    a.clear();
  }

  void fill(uint64_t i, uint64_t j, ValueT val);
  void insert(uint64_t i, uint64_t j, ValueT val);
  void remove(uint64_t i, uint64_t j, ValueT val);
  void defrag();

  void save(std::string filename, int file_format = Utils::BIN);
  void load(std::string filename);

  size_t size();

  bool empty();

  void clear();

  template <typename vecValueT>
  Vector<vecValueT> &operator*(Vector<vecValueT> &rvec);
  SparseMatrixCSR<ValueT> &operator*(ValueT rval);

  const std::string matrix_typename{"SparseMatrixDCSR"};

  using ValueType = ValueT;

private:
  std::vector<ValueT> a;
}; // class SparseMatrixDCSR

/*! Filling the sparse matrix with sorted entries by
 * ascending lines id
 */
template <typename ValueT>
void SparseMatrixDCSR<ValueT>::fill(uint64_t i, uint64_t j, ValueT val)
{
  SparseBMatrixDCSR::fill(i, j);
  a[ia[2 * pitch * i + 1]] = val;
}

/*! Inserting element in a DCSR matrix */
template <typename ValueT>
void SparseMatrixDCSR<ValueT>::insert(uint64_t i, uint64_t j, ValueT val)
{
}

/*! Remove element of the DCSR matrix */
template <typename ValueT>
void SparseMatrixDCSR<ValueT>::remove(uint64_t i, uint64_t j, ValueT val)
{
}

/*! Remove element of the DCSR matrix */
template <typename ValueT>
void SparseMatrixDCSR<ValueT>::defrag()
{
}

template <typename ValueT>
void SparseMatrixDCSR<ValueT>::save(std::string name, int fileformat)
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

  /* DCSR matrix properties */
  matfp.write(reinterpret_cast<const char *>(&pitch), sizeof(uint64_t));
  matfp.write(reinterpret_cast<const char *>(&alpha), sizeof(uint64_t));

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
void SparseMatrixDCSR<ValueT>::load(std::string name)
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

  /* DCSR matrix properties */
  matfp.read(reinterpret_cast<char *>(&pitch), sizeof(uint64_t));
  matfp.read(reinterpret_cast<char *>(&alpha), sizeof(uint64_t));

  if ((nnz + alpha + 2 * m * pitch) * sizeof(uint64_t)
      + (nnz + alpha) * sizeof(ValueT) < props->ram_limit)
  {
    a.resize(nnz + alpha, 0.);
    ia.resize(2 * pitch * m + 1, 0);
    ja.resize(nnz + alpha, 0);
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
size_t SparseMatrixDCSR<ValueT>::size()
{
  return (nnz + alpha + 2 * m * pitch) * sizeof(uint64_t)
         + (nnz + alpha) * sizeof(ValueT);
}


template <typename ValueT>
void SparseMatrixDCSR<ValueT>::clear()
{
  SparseBMatrixDCSR::clear();
  a.clear();
}

template <typename ValueT>
template <typename vecValueT>
Vector<vecValueT> &SparseMatrixDCSR<ValueT>::operator*(Vector<vecValueT> &rvec)
{
}

template <typename ValueT>
SparseMatrixDCSR<ValueT> &SparseMatrixDCSR<ValueT>::operator*(ValueT rval)
{
  for (uint64_t i = 0; i < nnz; i++)
    a[i] *= rval;

  return *this;
}

} // namespace graphee

#endif // GRAPHEE_SPARSE_MATRIX_DCSR_HPP__
