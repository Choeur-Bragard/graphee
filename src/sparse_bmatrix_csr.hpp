#ifndef GRAPHEE_SPARSE_BMATRIX_CSR_HPP__
#define GRAPHEE_SPARSE_BMATRIX_CSR_HPP__

#include <cstring>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include "snappy/snappy.h"

#include "properties.hpp"
#include "utils.hpp"
#include "vector.hpp"

namespace graphee {

/*! \brief Sparse boolean matrix in CSR format
 *         Commonly named Compressed Sparse Row
 *
 * The matrix only saves position of non-zero elements.
 * This is straightforward in non-weighted-edge graphs.
 */

class SparseBMatrixCSR {
public:
  SparseBMatrixCSR(Properties *properties)
      : props(properties), m(0), n(0), nnz(0) {}

  SparseBMatrixCSR(Properties *properties, uint64_t nlines, uint64_t ncols,
                   uint64_t nonzero_elems)
      : props(properties), m(nlines), n(ncols), nnz(nonzero_elems) {
    if ((nnz + m + 1) * sizeof(uint64_t) < props->ram_limit) {
      ia.resize(m + 1, 0);
      ja.resize(nnz, 0);
    } else {
      print_error("Requested size is beyond \'ram_limit\'");
      exit(-1);
    }
  }

  SparseBMatrixCSR(SparseBMatrixCSR &&mat)
      : props(mat.props), m(mat.m), n(mat.n), nnz(mat.nnz),
        ia(std::move(mat.ia)), ja(std::move(mat.ja)) {
    mat.props = nullptr;
    mat.m = 0;
    mat.n = 0;
    mat.nnz = 0;
  }

  SparseBMatrixCSR &operator=(SparseBMatrixCSR &&rmat);

  ~SparseBMatrixCSR() {
    ia.clear();
    ja.clear();
  }

  void fill(uint64_t i, uint64_t j);
  void insert(uint64_t i, uint64_t j);
  void remove(uint64_t i, uint64_t j);

  void save(std::string filename, int file_format = Utils::BIN);
  void load(std::string filename);

  size_t size();
  bool verify();

  bool empty();

  void clear();

  template <typename vecValueT>
  Vector<vecValueT> operator*(const Vector<vecValueT> &rvec);

  Vector<float> columns_sum();

  const std::string matrix_typename{"SparseBMatrixCSR"};

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
}; // class SparseBMatrixCSR

/*! Filling the sparse matrix with sorted entries by
 * ascending lines id
 */
void SparseBMatrixCSR::fill(uint64_t i, uint64_t j) {
  for (uint64_t l = fill_id + 1; l <= i; l++)
    ia[l + 1] = ia[l];

  ja[ia[i + 1]] = j;

  ia[i + 1]++;

  fill_id = i;
}

/*! Inserting element in a CSR matrix */
void SparseBMatrixCSR::insert(uint64_t i, uint64_t j) {}

/*! Remove element of the CSR matrix*/
void SparseBMatrixCSR::remove(uint64_t i, uint64_t j) {}

void SparseBMatrixCSR::save(std::string name, int fileformat) {
  std::ofstream matfp(name, std::ios_base::binary);

  size_t matrix_typename_size = matrix_typename.size();

  /* Save explicitly matrix properties */
  matfp.write(reinterpret_cast<const char *>(&matrix_typename_size),
              sizeof(size_t));
  matfp.write(reinterpret_cast<const char *>(matrix_typename.c_str()),
              matrix_typename_size);

  /* Save fileformat {BIN, SNAPPY} */
  matfp.write(reinterpret_cast<const char *>(&fileformat), sizeof(int));

  /* Matrix dimension */
  matfp.write(reinterpret_cast<const char *>(&m), sizeof(uint64_t));
  matfp.write(reinterpret_cast<const char *>(&nnz), sizeof(uint64_t));

  if (fileformat == Utils::BIN) {
    matfp.write(reinterpret_cast<const char *>(ia.data()),
                ia.size() * sizeof(uint64_t));
    matfp.write(reinterpret_cast<const char *>(ja.data()),
                ja.size() * sizeof(uint64_t));
  } else if (fileformat == Utils::SNAPPY) {
    size_t ia_snappy_size =
        snappy::MaxCompressedLength64(ia.size() * sizeof(uint64_t));
    char *ia_snappy = new char[ia_snappy_size];
    snappy::RawCompress64(reinterpret_cast<char *>(ia.data()),
                          ia.size() * sizeof(uint64_t), ia_snappy,
                          &ia_snappy_size);

    matfp.write(reinterpret_cast<const char *>(&ia_snappy_size),
                sizeof(size_t));
    matfp.write(reinterpret_cast<const char *>(ia_snappy), ia_snappy_size);
    delete[] ia_snappy;
    if (ja.size() != 0) {
      size_t ja_snappy_size =
          snappy::MaxCompressedLength64(ja.size() * sizeof(uint64_t));
      char *ja_snappy = new char[ja_snappy_size];
      snappy::RawCompress64(reinterpret_cast<char *>(ja.data()),
                            ja.size() * sizeof(uint64_t), ja_snappy,
                            &ja_snappy_size);

      matfp.write(reinterpret_cast<const char *>(&ja_snappy_size),
                  sizeof(size_t));
      matfp.write(reinterpret_cast<const char *>(ja_snappy), ja_snappy_size);
      delete[] ja_snappy;
    }
  }

  matfp.close();
}

void SparseBMatrixCSR::load(std::string name) {
  std::ifstream matfp(name, std::ios_base::binary);

  /* Save explicitly matrix properties */
  size_t matrix_typename_size;
  matfp.read(reinterpret_cast<char *>(&matrix_typename_size), sizeof(size_t));

  char read_matrix_typename[matrix_typename_size + 1];
  matfp.read(reinterpret_cast<char *>(read_matrix_typename),
             matrix_typename_size);
  read_matrix_typename[matrix_typename_size] = '\0';

  if (std::strcmp(read_matrix_typename, matrix_typename.c_str()) != 0) {
    std::ostringstream oss;
    oss << "Wrong matrix format, found \'" << read_matrix_typename
        << "\' while expecting \'" << matrix_typename << "\'";
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

  if ((nnz + m + 1) * sizeof(uint64_t) < props->ram_limit) {
    ia.resize(m + 1, 0);
    ja.resize(nnz, 0);
  } else {
    print_error("Requested size is beyond \'ram_limit\'");
    matfp.close();
    exit(-1);
  }

  if (fileformat == Utils::BIN) {
    matfp.read(reinterpret_cast<char *>(ia.data()),
               ia.size() * sizeof(uint64_t));
    matfp.read(reinterpret_cast<char *>(ja.data()),
               ja.size() * sizeof(uint64_t));
  } else if (fileformat == Utils::SNAPPY) {
    bool uncomp_succeed;

    size_t ia_snappy_size;
    matfp.read(reinterpret_cast<char *>(&ia_snappy_size), sizeof(size_t));

    char *ia_snappy = new char[ia_snappy_size];
    matfp.read(reinterpret_cast<char *>(ia_snappy), ia_snappy_size);

    uncomp_succeed = snappy::RawUncompress64(
        ia_snappy, ia_snappy_size, reinterpret_cast<char *>(ia.data()));
    delete[] ia_snappy;

    if (!uncomp_succeed) {
      print_error("SNAPPY uncompression of IA failed");

      matfp.close();
      exit(-1);
    }

    if (nnz != 0) {
      size_t ja_snappy_size;
      matfp.read(reinterpret_cast<char *>(&ja_snappy_size), sizeof(size_t));

      char *ja_snappy = new char[ja_snappy_size];
      matfp.read(reinterpret_cast<char *>(ja_snappy), ja_snappy_size);

      uncomp_succeed = snappy::RawUncompress64(
          ja_snappy, ja_snappy_size, reinterpret_cast<char *>(ja.data()));
      delete[] ja_snappy;

      if (!uncomp_succeed) {
        print_error("SNAPPY uncompression of JA failed");

        matfp.close();
        exit(-1);
      }
    }
  }

  fill_id = m - 1;

  matfp.close();
}

size_t SparseBMatrixCSR::size() { return (nnz + m + 1) * sizeof(uint64_t); }

bool SparseBMatrixCSR::verify() {
  if (fill_id < m - 1) {
    for (uint64_t l = fill_id + 1; l < m; l++) {
      ia[l + 1] = ia[l];
    }
    fill_id = m - 1;
  }

  if (nnz == ia[m]) {
    return true;
  } else {
    std::ostringstream oss;
    oss << "NNZ = " << nnz << " IA[M] = " << ia[m];
    print_warning(oss.str());
    return false;
  }
}

void SparseBMatrixCSR::clear() {
  ia.clear();
  ja.clear();

  m = 0;
  n = 0;
  nnz = 0;
}

SparseBMatrixCSR &SparseBMatrixCSR::operator=(SparseBMatrixCSR &&rmat) {
  std::swap(ia, rmat.ia);
  std::swap(ja, rmat.ja);
  std::swap(m, rmat.m);
  std::swap(n, rmat.n);
  std::swap(nnz, rmat.nnz);

  return *this;
}

// template <typename vecValueT=float>
Vector<float> SparseBMatrixCSR::columns_sum() {

  Vector<float> res(props, m, 0.);

#pragma omp parallel for num_threads(props->nthreads)
  for (uint64_t i = 0; i < ja.size(); i++) {
    res[ja[i]] += 1;
  }

  return res;
}

template <typename vecValueT>
Vector<vecValueT> SparseBMatrixCSR::operator*(const Vector<vecValueT> &rvec) {
  if (n != rvec.get_lines()) {
    std::ostringstream oss;
    oss << "Error SpBMat[" << m << "x" << n << "] with Vec[" << rvec.get_lines()
        << "]";
    print_error(oss.str());
    exit(-1);
  }

  Vector<vecValueT> res(props, m, 0.);

#pragma omp parallel for num_threads(props->nthreads)
  for (uint64_t i = 0; i < m; i++) {
    for (uint64_t ja_idx = ia[i]; ja_idx < ia[i + 1]; ja_idx++) {
      res[i] += rvec[ja[ja_idx]];
    }
  }

  return res;
}

uint64_t SparseBMatrixCSR::get_lines() { return m; }

uint64_t SparseBMatrixCSR::get_columns() { return n; }

uint64_t SparseBMatrixCSR::get_nonzeros() { return nnz; }

} // namespace graphee

#endif // GRAPHEE_SPARSE_BMATRIX_CSR_HPP__
