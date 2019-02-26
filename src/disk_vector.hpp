#ifndef GRAPHEE_DISKVEC_HPP__
#define GRAPHEE_DISKVEC_HPP__

#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include <cstdio>

#include "disk_sparse_matrix.hpp"
#include "vector.hpp"

/* \brief Vector saved by slices within disk
 *
 * Modify the format in order to template a vector
 * standart class instead of a value type.
 * In such cases either sparse or dense vector could be written,
 * on the disk.
 *
 */

namespace graphee {

template <typename VectorT> class DiskVector {
public:
  DiskVector() : props(nullptr) {}

  DiskVector(Properties *properties) : props(properties) {}

  DiskVector(Properties *properties, std::string vector_name,
             typename VectorT::ValueType init_val = 0)
      : props(properties), name(vector_name), m(properties->nvertices) {
    if (props->window * sizeof(typename VectorT::ValueType) >
        props->ram_limit) {
      print_error("The \'graphee::vector\' size exceeds the \'ram_limit\'");
      exit(-1);
    }

    VectorT tmp(props, props->window, init_val);
    for (uint64_t slice_id = 0; slice_id < props->nslices; slice_id++) {
      tmp.save(get_slice_filename(slice_id));
    }
  }

  DiskVector(DiskVector<VectorT> &&vec)
      : props(vec.props), name(vec.name), m(vec.props->nvertices) {
    vec.props = nullptr;
    vec.name = "";
    vec.m = 0;
  }

  ~DiskVector() {}

  uint64_t get_n_slices();
  VectorT &get_slice(uint64_t slice_id);

  void swap(DiskVector<VectorT> &rvec);

  template <typename DiskMatrixT> void dmat_columns_sum(DiskMatrixT &dmat);

  uint64_t countZeros();

  DiskVector<VectorT> &operator=(DiskVector &&vec) {
    std::swap(props, vec.props);
    std::swap(name, vec.name);
    std::swap(m, vec.m);
  }

  template <typename DiskMatrixT>
  void dmat_prod_dvec(typename VectorT::ValueType a, DiskMatrixT &dmat,
                      DiskVector<VectorT> &dvec);

  template <typename DiskMatrixT>
  void dmat_prod_dvec_over_dvec(typename VectorT::ValueType a,
                                DiskMatrixT &dmat, DiskVector<VectorT> &ldvec,
                                DiskVector<VectorT> &rdvec);

  /* Divide each term of this by the term of ldvec with the same index and
   * accumulate values of this when the division is undefined : If ldvec[i] != 0
   * then this[i]/=ldvec[i] If ldvec[i] == 0 then aggregator+=this[i] ;
   * this[i]=0;
   */
  void divide_and_sum_Nan(DiskVector<VectorT> &ldvec, float &aggregator);

  DiskVector<VectorT> &operator+=(typename VectorT::ValueType val);

  using VectorType = VectorT;

  uint64_t m;

private:
  Properties *props;
  std::string name;

  std::string get_slice_filename(uint64_t slice_id);
}; // class DiskVector

template <typename VectorT> uint64_t DiskVector<VectorT>::get_n_slices() {
  return props->nslices;
}

template <typename VectorT>
VectorT &DiskVector<VectorT>::get_slice(uint64_t slice_id) {
  std::ostringstream oss;
  oss << "Load slice [" << slice_id << "] of vector \'" << name << "\'";
  print_log(oss.str());

  VectorT &res = *new VectorT(props);
  res.load(get_slice_filename(slice_id));

  return res;
}

template <typename VectorT>
std::string DiskVector<VectorT>::get_slice_filename(uint64_t slice_id) {
  std::ostringstream slicename;
  slicename << props->name << "_" << name << "_dvecslc_" << slice_id << ".gpe";
  return slicename.str();
}

template <typename VectorT>
void DiskVector<VectorT>::swap(DiskVector<VectorT> &vec) {
  if (props->nvertices != vec.props->nvertices) {
    std::ostringstream oss;
    oss << "Could not swap vector \'" << name << "\' and \'";
    oss << vec.name << "\' because dimensions are not equal";
    print_error(oss.str());
  }

  int succed;
  std::ostringstream tmpname;
  tmpname << name << "_swap_file.gpe";
  for (uint64_t slice_id = 0; slice_id < props->nslices; slice_id++) {
    succed =
        rename(get_slice_filename(slice_id).c_str(), tmpname.str().c_str());

    if (succed != 0) {
      std::ostringstream oss;
      oss << "Could not swap vector \'" << get_slice_filename(slice_id)
          << "\' to \'";
      oss << tmpname.str() << "\'";
      print_error(oss.str());
    }

    succed = rename(vec.get_slice_filename(slice_id).c_str(),
                    get_slice_filename(slice_id).c_str());

    if (succed != 0) {
      std::ostringstream oss;
      oss << "Could not swap vector \'" << vec.get_slice_filename(slice_id)
          << "\' to \'";
      oss << get_slice_filename(slice_id) << "\'";
      print_error(oss.str());
    }

    succed =
        rename(tmpname.str().c_str(), vec.get_slice_filename(slice_id).c_str());

    if (succed != 0) {
      std::ostringstream oss;
      oss << "Could not swap vector \'" << tmpname.str() << "\' to \'";
      oss << vec.get_slice_filename(slice_id) << "\'";
      print_error(oss.str());
    }
  }
}

template <typename VectorT>
template <typename DiskMatrixT>
void DiskVector<VectorT>::dmat_columns_sum(DiskMatrixT &dmat) {
  if (dmat.n != m) {
    std::ostringstream oss;
    oss << "Wrong dimensions in matrix sum columns";
    print_error(oss.str());
    exit(-1);
  }
#pragma omp parallel for
  for (uint64_t col = 0; col < props->nslices; col++) {
    VectorT res(std::move(this->get_slice(col)));
    for (uint64_t line = 0; line < props->nslices; line++) {
      typename DiskMatrixT::MatrixType smat = std::move(dmat.get_block(line, col));
      
      res += smat.columns_sum(); 
    }
    res.save(this->get_slice_filename(col));
  }
}

template <typename VectorT> uint64_t DiskVector<VectorT>::countZeros() {
  uint64_t res = 0;
#pragma omp parallel for reduction(+ : res)
  for (uint64_t slice = 0; slice < props->nslices; slice++) {
    VectorT vectSlice(std::move(this->get_slice(slice)));
    for (uint64_t i = 0; i < props->window; i++) {
      if (!vectSlice[i]) {
        res += 1;
      }
    }
  }
  return res;
}

template <typename VectorT>
DiskVector<VectorT> &DiskVector<VectorT>::
operator+=(typename VectorT::ValueType val) {
  VectorT vec(props);
  for (uint64_t slice_id = 0; slice_id < props->nslices; slice_id++) {
    vec.load(get_slice_filename(slice_id));
    vec += val;
    vec.save(get_slice_filename(slice_id));
  }

  return (*this);
}

template <typename VectorT>
template <typename DiskMatrixT>
void DiskVector<VectorT>::dmat_prod_dvec(typename VectorT::ValueType a,
                                         DiskMatrixT &dmat,
                                         DiskVector<VectorT> &dvec) {
  if (dmat.n != dvec.m) {
    std::ostringstream oss;
    oss << "Wrong dimensions in Matrix-Vector product";
    print_error(oss.str());
    exit(-1);
  }

#pragma omp parallel for
  for (uint64_t line = 0; line < props->nslices; line++) {
    VectorT res(std::move(this->get_slice(line)));
    for (uint64_t col = 0; col < props->nslices; col++) {
      typename DiskMatrixT::MatrixType smat =
          std::move(dmat.get_block(line, col));
      VectorT rvec(std::move(dvec.get_slice(col)));

      rvec *= a;
      res += smat * rvec;
    }
    res.save(this->get_slice_filename(line));
  }
}

template <typename VectorT>
void DiskVector<VectorT>::divide_and_sum_Nan(DiskVector<VectorT> &ldvec,
                                             float &aggregator) {
  if (props->nslices != ldvec.props->nslices) {
    std::ostringstream oss;
    oss << "Wrong dimensions in Hadamard division, both vector must have the "
           "same size";
    print_error(oss.str());
    exit(-1);
  }

// float aggregator=0;
#pragma omp parallel for reduction(+ : aggregator)
  for (uint64_t slice = 0; slice < props->nslices; slice++) {
    VectorT res(std::move(this->get_slice(slice)));
    VectorT lvec(std::move(ldvec.get_slice(slice)));
    res.divide_and_sum_Nan(lvec, aggregator);
    res.save(this->get_slice_filename(slice));
  }
}

template <typename VectorT>
template <typename DiskMatrixT>
void DiskVector<VectorT>::dmat_prod_dvec_over_dvec(
    typename VectorT::ValueType a, DiskMatrixT &dmat,
    DiskVector<VectorT> &ldvec, DiskVector<VectorT> &rdvec) {
  if (dmat.n != ldvec.m && ldvec.m != rdvec.m) {
    std::ostringstream oss;
    oss << "Wrong dimensions in Matrix-Vector product";
    print_error(oss.str());
    exit(-1);
  }

#pragma omp parallel for
  for (uint64_t line = 0; line < props->nslices; line++) {
    VectorT res(std::move(this->get_slice(line)));
    for (uint64_t col = 0; col < props->nslices; col++) {
      typename DiskMatrixT::MatrixType smat =
          std::move(dmat.get_block(line, col));
      VectorT lvec(std::move(ldvec.get_slice(col)));
      VectorT rvec(std::move(rdvec.get_slice(col)));

      lvec /= rvec;
      lvec *= a;
      res += smat * lvec;
    }
    res.save(this->get_slice_filename(line));
  }
}

} // namespace graphee

#endif // GRAPHEE_DISKVEC_H__
