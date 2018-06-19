#ifndef GRAPHEE_SPARSE_BMATRIX_DCSR_HPP__
#define GRAPHEE_SPARSE_BMATRIX_DCSR_HPP__

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cstring>
#include <vector>

#include "snappy/snappy.h"

#include "utils.hpp"
#include "properties.hpp"
#include "vector.hpp"

namespace graphee
{

/*! \brief Sparse boolean matrix in Dynamic CSR format
 *         Commonly named Dynamic Compressed Sparse Row:
 *         https://thomas.gilray.org/pdf/dynamic-csr.pdf
 *
 * The matrix only saves position of non-zero elements.
 * This is straightforward in non-weighted-edge graphs.
 */

class SparseBMatrixCSR
{
public:
  SparseBMatrixDCSR(Properties *properties) :
    props(properties), m(0), n(0), nnz(0), pitch(0), alpha(0) {}

  SparseBMatrixDCSR(Properties *properties, uint64_t nlines, uint64_t ncols,
                    uint64_t nonzero_elems, uint64_t pitch, uint64_t alpha) :
    props(properties), m(nlines), n(ncols), nnz(nonzero_elems), pitch(pitch),
    alpha(alpha)
  {
    if ((nnz + alpha + 2 * m * pitch) * sizeof(uint64_t) < props->ram_limit)
    {
      ia.resize(2 * m * pitch, 0);
      ja.resize(nnz + alpha, 0);
    }
    else
    {
      print_error("Requested size is beyond \'ram_limit\'");
      exit(-1);
    }
  }

  SparseBMatrixCSR(SparseBMatrixCSR &&mat) : props(mat.props), m(mat.m), n(mat.n),
    nnz(mat.nnz), ia(std::move(mat.ia)), ja(std::move(mat.ja)), pitch(mat.pitch),
    alpha(mat.alpha)
  {
    mat.props = nullptr;
    mat.m = 0;
    mat.n = 0;
    mat.nnz = 0;
    mat.pitch = 0;
    mat.alpha = 0;
  }

  SparseBMatrixCSR &operator=(SparseBMatrixCSR &&rmat);

  ~SparseBMatrixCSR()
  {
    ia.clear();
    ja.clear();
  }

  void fill(uint64_t i, uint64_t j);
  void insert(uint64_t i, uint64_t j);
  void remove(uint64_t i, uint64_t j);
  void defrag();

  void save(std::string filename, int file_format = Utils::BIN);
  void load(std::string filename);

  size_t size();
  bool verify();

  bool empty();

  void clear();

  template <typename vecValueT>
  Vector<vecValueT> operator*(const Vector<vecValueT> &rvec);

  const std::string matrix_typename{"SparseBMatrixDCSR"};

  uint64_t get_lines();
  uint64_t get_columns();
  uint64_t get_nonzeros();

private:
  std::vector<uint64_t> ia;
  std::vector<uint64_t> ja;

  Properties *props;

  uint64_t m;
  uint64_t n;
  uint64_t nnz;
  uint64_t fill_id;

  uint64_t alpha;
  uint64_t pitch;
}; // class SparseBMatrixDCSR

/*! Filling the sparse matrix with sorted entries by
 * ascending lines id
 */
void SparseBMatrixDCSR::fill(uint64_t i, uint64_t j)
{
  for (uint64_t l = 2 * pitch * fill_id + 1; // first ID
       l < 2 * pitch * i + 1; // last ID
       l += 2 * pitch)
  {
    ia[l + 2 * pitch] = ia[l];
    ia[l + 2 * pitch - 1] = ia[l + 2 * pitch]
  }

  ja[ia[2 * pitch * i + 1]] = j;

  ia[2 * pitch * i + 1]++;

  fill_id = i;
}

/*! Inserting element in a DCSR matrix */
void SparseBMatrixDCSR::insert(uint64_t i, uint64_t j)
{
}

/*! Remove element of the DCSR matrix */
void SparseBMatrixDCSR::remove(uint64_t i, uint64_t j)
{
}

/*! Defragment the DCSR matrix */
void SparseBMatrixDCSR::defrag()
{
}

void SparseBMatrixDCSR::save(std::string name, int fileformat)
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

  /* DCSR properties */
  matfp.write(reinterpret_cast<const char *>(&pitch), sizeof(uint64_t));
  matfp.write(reinterpret_cast<const char *>(&alpha), sizeof(uint64_t));

  if (fileformat == Utils::BIN)
  {
    matfp.write(reinterpret_cast<const char *>(ia.data()), ia.size() * sizeof(uint64_t));
    matfp.write(reinterpret_cast<const char *>(ja.data()), ja.size() * sizeof(uint64_t));
  }
  else if (fileformat == Utils::SNAPPY)
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

void SparseBMatrixDCSR::load(std::string name)
{
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
  n = m;

  /* DCSR matrix properties */
  matfp.read(reinterpret_cast<char *>(&pitch), sizeof(uint64_t));
  matfp.read(reinterpret_cast<char *>(&alpha), sizeof(uint64_t));

  if ((nnz + alpha + 2 * m * pitch) * sizeof(uint64_t) < props->ram_limit)
  {
    ia.resize(2 * pitch * m, 0);
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
    matfp.read(reinterpret_cast<char *>(ia.data()), ia.size() * sizeof(uint64_t));
    matfp.read(reinterpret_cast<char *>(ja.data()), ja.size() * sizeof(uint64_t));
  }
  else if (fileformat == Utils::SNAPPY)
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

size_t SparseBMatrixDCSR::size()
{
  return (nnz + alpha + 2 * pitch * m + 1) * sizeof(uint64_t);
}

bool SparseBMatrixDCSR::verify()
{
  for (uint64_t l = 2 * pitch * fill_id + 1; // first ID
       l < 2 * pitch * (m - 1) + 1; // last ID
       l += 2 * pitch)
  {
    ia[l + 2 * pitch] = ia[l];
    ia[l + 2 * pitch - 1] = ia[l + 2 * pitch]
  }
  fill_id = m - 1;

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

void SparseBMatrixDCSR::clear()
{
  ia.clear();
  ja.clear();

  m = 0;
  n = 0;
  nnz = 0;
  pitch = 0;
  alpha = 0;
}

SparseBMatrixDCSR &SparseBMatrixDCSR::operator=(SparseBMatrixDCSR &&rmat)
{
  std::swap(ia, rmat.ia);
  std::swap(ja, rmat.ja);
  std::swap(m, rmat.m);
  std::swap(n, rmat.n);
  std::swap(pitch, rmat.pitch);
  std::swap(alpha, rmat.alpha);
  std::swap(nnz, rmat.nnz);

  return *this;
}

template <typename vecValueT>
Vector<vecValueT> SparseBMatrixDCSR::operator*(const Vector<vecValueT> &rvec)
{
}

uint64_t SparseBMatrixDCSR::get_lines()
{
  return m;
}

uint64_t SparseBMatrixDCSR::get_columns()
{
  return n;
}

uint64_t SparseBMatrixDCSR::get_nonzeros()
{
  return nnz;
}

} // namespace graphee

#endif // GRAPHEE_SPARSE_BMATRIX_DCSR_HPP__
