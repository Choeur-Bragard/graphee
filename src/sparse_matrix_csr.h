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
#include "sparse_bmatrix_csr.h"

namespace graphee
{

/*! \brief Sparse matrix in CSR format
 *         Commonly named Compressed Sparse Row
 *
 * The matrix can handle either `double`, `float` and `int` 
 * values.
 */

template <typename valueT>
class sparseMatrixCSR : public sparseBMatrixCSR
{
public:
  sparseMatrixCSR(properties *properties) : sparseBMatrixCSR(properties)
  {
    if (typeid(valueT) == typeid(bool))
    {
      print_error ("Do not use \'sparseMatrixCSR<bool>\', instead use \'sparseBMatrixCSR\'");
    }
  }

  sparseMatrixCSR(properties *properties, uint64_t nlines, uint64_t ncols,
                  uint64_t nonzero_elems, valueT init_val = 0.) :
                  sparseBMatrixCSR(properties, nlines, ncols, nonzero_elems)
  {
    if (typeid(valueT) == typeid(bool))
    {
      print_error ("Do not use \'sparseMatrixCSR<bool>\', instead use \'sparseBMatrixCSR\'");
    }

    if ((nnz + m + 1) * sizeof(uint64_t) + nnz * sizeof(valueT) < props->ram_limit)
    {
      a.resize(nnz, init_val);
    }
    else
    {
      print_error("Requested size is beyond \'ram_limit\'");
      exit(-1);
    }
  }

  ~sparseMatrixCSR() {}

  void fill(uint64_t i, uint64_t j, valueT val);
  void insert(uint64_t i, uint64_t j, valueT val);
  void remove(uint64_t i, uint64_t j, valueT val);

  void save(std::string filename, int file_format = utils::BIN);
  void load(std::string filename);

  size_t size();

  bool empty();

  void clear();

  template <typename vecValueT>
  vector<vecValueT> &operator*(vector<vecValueT> &rvec);
  sparseMatrixCSR<valueT> &operator*(valueT rval);

  const std::string matrixType{"sparseMatrixCSR"};

  using valueType = valueT;

private:
  std::vector<valueT> a;
}; // class sparseMatrixCSR

/*! Filling the sparse matrix with sorted entries by
 * ascending lines id
 */
template <typename valueT>
void sparseMatrixCSR<valueT>::fill(uint64_t i, uint64_t j, valueT val)
{
  sparseBMatrixCSR::fill(i, j);
  a[ia[i + 1]] = val;
}

/*! Inserting element in a CSR matrix */
template <typename valueT>
void sparseMatrixCSR<valueT>::insert(uint64_t i, uint64_t j, valueT val)
{
}

/*! Remove element of the CSR matrix*/
template <typename valueT>
void sparseMatrixCSR<valueT>::remove(uint64_t i, uint64_t j, valueT val)
{
}

template <typename valueT>
void sparseMatrixCSR<valueT>::save(std::string name, int fileformat)
{
  std::ofstream matfp(name, std::ios_base::binary);

  size_t matrixType_size = matrixType.size();

  /* Save explicitly matrix properties */
  matfp.write(reinterpret_cast<const char *>(&matrixType_size), sizeof(size_t));
  matfp.write(reinterpret_cast<const char *>(matrixType.c_str()), matrixType_size);

  /* Save fileformat {BIN, SNAPPY} */
  matfp.write(reinterpret_cast<const char *>(&fileformat), sizeof(int));

  /* Matrix dimension */
  matfp.write(reinterpret_cast<const char *>(&m), sizeof(uint64_t));
  matfp.write(reinterpret_cast<const char *>(&nnz), sizeof(uint64_t));

  if (fileformat == utils::BIN)
  {
    matfp.write(reinterpret_cast<const char *>(a.data()), a.size() * sizeof(valueT));
    matfp.write(reinterpret_cast<const char *>(ia.data()), ia.size() * sizeof(uint64_t));
    matfp.write(reinterpret_cast<const char *>(ja.data()), ja.size() * sizeof(uint64_t));
  }
  else if (fileformat == utils::SNAPPY)
  {
    size_t a_snappy_size = snappy::MaxCompressedLength64(a.size() * sizeof(valueT));
    char *a_snappy = new char[a_snappy_size];
    snappy::RawCompress64(reinterpret_cast<char *>(a.data()), a.size() * sizeof(valueT), a_snappy, &a_snappy_size);

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

template <typename valueT>
void sparseMatrixCSR<valueT>::load(std::string name)
{
  this->clear();
  std::ifstream matfp(name, std::ios_base::binary);

  /* Save explicitly matrix properties */
  size_t matrixType_size;
  matfp.read(reinterpret_cast<char *>(&matrixType_size), sizeof(size_t));

  char read_matrixType[matrixType_size + 1];
  matfp.read(reinterpret_cast<char *>(read_matrixType), matrixType_size);
  read_matrixType[matrixType_size] = '\0';

  if (std::strcmp(read_matrixType, matrixType.c_str()) != 0)
  {
    std::ostringstream oss;
    oss << "Wrong matrix format, found \'" << read_matrixType << "\' while expecting \'"
        << matrixType << "\'";
    print_error(oss.str());
    exit(-1);
  }

  /* Read fileformat {BIN, SNAPPY} */
  int fileformat;
  matfp.read(reinterpret_cast<char *>(&fileformat), sizeof(int));

  /* Matrix dimension */
  matfp.read(reinterpret_cast<char *>(&m), sizeof(uint64_t));
  matfp.read(reinterpret_cast<char *>(&nnz), sizeof(uint64_t));

  if ((nnz + m + 1) * sizeof(uint64_t) + nnz * sizeof(valueT) < props->ram_limit)
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

  if (fileformat == utils::BIN)
  {
    matfp.read(reinterpret_cast<char *>(a.data()), a.size() * sizeof(valueT));
    matfp.read(reinterpret_cast<char *>(ia.data()), ia.size() * sizeof(uint64_t));
    matfp.read(reinterpret_cast<char *>(ja.data()), ja.size() * sizeof(uint64_t));
  }
  else if (fileformat == utils::SNAPPY)
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

template <typename valueT>
size_t sparseMatrixCSR<valueT>::size()
{
  return (nnz + m + 1) * sizeof(uint64_t) + nnz * sizeof(valueT);
}


template <typename valueT>
void sparseMatrixCSR<valueT>::clear()
{
  sparseBMatrixCSR::clear();
  a.clear();
}

template <typename valueT>
template <typename vecValueT>
vector<vecValueT> &sparseMatrixCSR<valueT>::operator*(vector<vecValueT> &rvec)
{
  if (n != rvec.get_lines())
  {
    std::ostringstream oss;
    oss << "Error SpMat[" << m << "x" << n << "] with Vec[" << rvec.get_lines() << "]";
    print_error(oss.str());
    exit(-1);
  }

  vector<vecValueT> res(rvec.get_properties(), m, 0.);

  for (uint64_t i = 0; i < m; i++)
  {
    for (uint64_t ja_idx = ia[i]; ja_idx < ia[i + 1]; ja_idx++)
    {
      res[i] += a[ja[ja_idx]] * rvec[ja[ja_idx]];
    }
  }
}

template <typename valueT>
sparseMatrixCSR<valueT> &sparseMatrixCSR<valueT>::operator*(valueT rval)
{
  for (uint64_t i = 0; i < nnz; i++)
    a[i] *= rval;

  return *this;
}

} // namespace graphee

#endif // GRAPHEE_SPARSE_MATRIX_CSR_H__
