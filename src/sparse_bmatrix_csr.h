#ifndef GRAPHEE_SPARSE_BMATRIX_CSR_H__
#define GRAPHEE_SPARSE_BMATRIX_CSR_H__

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

namespace graphee
{

/*! \brief Sparse boolean matrix in CSR format
 *         Commonly named Compressed Sparse Row
 *
 * The matrix only saves position of non-zero elements.
 * This is straightforward in non-weighted-edge graphs.
 */

class sparseBMatrixCSR
{
public:
  sparseBMatrixCSR(properties *properties) : props(properties), m(0), n(0), nnz(0)
  {
    std::cout << "SpBMat empty" << std::endl;
  }

  sparseBMatrixCSR(properties *properties, uint64_t nlines, uint64_t ncols,
                   uint64_t nonzero_elems) : props(properties), m(nlines),
    n(ncols), nnz(nonzero_elems)
  {
    std::cout << "SpBMat full" << std::endl;
    if ((nnz + m + 1) * sizeof(uint64_t) < props->ram_limit)
    {
      ia.resize(m + 1, 0);
      ja.resize(nnz, 0);
    }
    else
    {
      print_error("Requested size is beyond \'ram_limit\'");
      exit(-1);
    }
  }

  ~sparseBMatrixCSR()
  {
    std::cout << "SpBMat destructor" << std::endl;
    ia.clear();
    ja.clear();
  }

  void fill(uint64_t i, uint64_t j);
  void insert(uint64_t i, uint64_t j);
  void remove(uint64_t i, uint64_t j);

  void save(std::string filename, int file_format = utils::BIN);
  void load(std::string filename);

  size_t size();
  bool verify();

  bool empty();

  void clear();

  template <typename vecValueT>
  vector<vecValueT> &operator*(vector<vecValueT> &rvec);

  const std::string matrixType{"sparseBMatrixCSR"};

  uint64_t get_lines();
  uint64_t get_columns();
  uint64_t get_nonzeros();

private:
  std::vector<uint64_t> ia;
  std::vector<uint64_t> ja;

  properties *props;

  uint64_t m;
  uint64_t n;
  uint64_t nnz;
  uint64_t fill_id;
}; // class sparseBMatrixCSR

/*! Filling the sparse matrix with sorted entries by
 * ascending lines id
 */
void sparseBMatrixCSR::fill(uint64_t i, uint64_t j)
{
  for (uint64_t l = fill_id + 1; l <= i; l++)
    ia[l + 1] = ia[l];

  ja[ia[i + 1]] = j;

  ia[i + 1]++;

  fill_id = i;
}

/*! Inserting element in a CSR matrix */
void sparseBMatrixCSR::insert(uint64_t i, uint64_t j)
{
}

/*! Remove element of the CSR matrix*/
void sparseBMatrixCSR::remove(uint64_t i, uint64_t j)
{
}

void sparseBMatrixCSR::save(std::string name, int fileformat)
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
    matfp.write(reinterpret_cast<const char *>(ia.data()), ia.size() * sizeof(uint64_t));
    matfp.write(reinterpret_cast<const char *>(ja.data()), ja.size() * sizeof(uint64_t));
  }
  else if (fileformat == utils::SNAPPY)
  {
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

void sparseBMatrixCSR::load(std::string name)
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
  n = m;

  if ((m + 1) * sizeof(uint64_t) + nnz * sizeof(uint64_t) < props->ram_limit)
  {
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
    matfp.read(reinterpret_cast<char *>(ia.data()), ia.size() * sizeof(uint64_t));
    matfp.read(reinterpret_cast<char *>(ja.data()), ja.size() * sizeof(uint64_t));
  }
  else if (fileformat == utils::SNAPPY)
  {
    bool uncomp_succeed;

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

size_t sparseBMatrixCSR::size()
{
  return (nnz + m + 1) * sizeof(uint64_t);
}

bool sparseBMatrixCSR::verify()
{
  if (fill_id < m - 1)
  {
    for (uint64_t l = fill_id + 1; l < m; l++)
    {
      ia[l + 1] = ia[l];
    }
    fill_id = m - 1;
  }

  if (nnz == ia[m])
  {
    return true;
  }
  else
  {
    std::ostringstream oss;
    oss << "NNZ = " << nnz << " IA[M] = " << ia[m];
    print_warning(oss.str());
    return false;
  }
}

void sparseBMatrixCSR::clear()
{
  ia.clear();
  ja.clear();

  m = 0;
  n = 0;
  nnz = 0;
}

template <typename vecValueT>
vector<vecValueT> &sparseBMatrixCSR::operator*(vector<vecValueT> &rvec)
{
  if (n != rvec.get_lines())
  {
    std::ostringstream oss;
    oss << "Error SpBMat[" << m << "x" << n << "] with Vec[" << rvec.get_lines() << "]";
    print_error(oss.str());
    exit(-1);
  }

  vector<vecValueT> &res = *new vector<vecValueT>(rvec.get_properties(), m, 0.);

  for (uint64_t i = 0; i < m; i++)
  {
    for (uint64_t ja_idx = ia[i]; ja_idx < ia[i + 1]; ja_idx++)
    {
      res[i] += rvec[ja[ja_idx]];
    }
  }

  return res;
}

uint64_t sparseBMatrixCSR::get_lines()
{
  return m;
}

uint64_t sparseBMatrixCSR::get_columns()
{
  return n;
}

uint64_t sparseBMatrixCSR::get_nonzeros()
{
  return nnz;
}

} // namespace graphee

#endif // GRAPHEE_SPARSE_BMATRIX_CSR_H__
