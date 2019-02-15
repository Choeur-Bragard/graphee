#ifndef GRAPHEE_VECTOR_HPP__
#define GRAPHEE_VECTOR_HPP__

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>

#include "utils.hpp"
#include "properties.hpp"

#include "sparse_bmatrix_csr.hpp"
#include "sparse_matrix_csr.hpp"

namespace graphee
{

/*! \brief Graphee Vector class
 *         Includes load/dump of data on disk
 *
 * It contains overloaded operators for operations
 * with `sparseMatrix*` classes
 */

template <typename ValueT>
class Vector : public std::vector<ValueT>
{
public:
  Vector(Properties *properties) :
    std::vector<ValueT>(), props(properties) {}

  Vector(Properties *properties, size_t m, ValueT init_value = 0) :
    std::vector<ValueT>(m, init_value), props(properties) {}

  Vector(Vector &&vec) : std::vector<ValueT>(std::move(vec)), props(vec.props) {}

  ~Vector() {}

  void save(std::string name, int fileformat = Utils::BIN);
  void load(std::string name);

  Vector<ValueT> &operator+=(const Vector<ValueT>& rvec);
  Vector<ValueT> &operator+=(ValueT val);
  Vector<ValueT> &operator*=(ValueT val);

  Vector<ValueT> &operator=(Vector<ValueT> &&rvec);

  Vector<ValueT> &operator/=(Vector<ValueT> &rvec);

  Vector<ValueT> &divide_and_sum_Nan(Vector<ValueT> &rvec, ValueT& aggs);

  const std::string vector_typename{"Vector"};

  uint64_t get_lines() const;

  using ValueType = ValueT;

  Properties* get_properties();

private:
  Properties *props;

}; // class graphee::Vector

template <class ValueT>
void Vector<ValueT>::save(std::string name, int fileformat)
{
  std::ofstream vecfp(name, std::ios_base::binary);

  size_t vector_typename_size = vector_typename.size();

  /* Save explicitly Vector properties */
  vecfp.write(reinterpret_cast<const char *>(&vector_typename_size), sizeof(size_t));
  vecfp.write(reinterpret_cast<const char *>(vector_typename.c_str()), vector_typename_size);

  /* Save fileformat {BIN, SNAPPY} */
  vecfp.write(reinterpret_cast<const char *>(&fileformat), sizeof(int));

  /* Vector dimension */
  size_t vec_size = this->size();
  vecfp.write(reinterpret_cast<const char *>(&vec_size), sizeof(size_t));

  if (fileformat == Utils::BIN)
  {
    vecfp.write(reinterpret_cast<const char *>(this->data()), this->size() * sizeof(ValueT));
  }
  else if (fileformat == Utils::SNAPPY)
  {
    size_t vec_snappy_size = snappy::MaxCompressedLength64(this->size() * sizeof(ValueT));
    char *vec_snappy = new char[vec_snappy_size];
    snappy::RawCompress64(reinterpret_cast<char *>(this->data()), this->size() * sizeof(ValueT), vec_snappy, &vec_snappy_size);

    vecfp.write(reinterpret_cast<const char *>(&vec_snappy_size), sizeof(size_t));
    vecfp.write(reinterpret_cast<const char *>(vec_snappy), vec_snappy_size);
    delete[] vec_snappy;
  }

  vecfp.close();
}

template <class ValueT>
void Vector<ValueT>::load(std::string name)
{
  this->clear();
  std::ifstream vecfp(name, std::ios_base::binary);

  /* Save explicitly Vector properties */
  size_t vector_typename_size;
  vecfp.read(reinterpret_cast<char *>(&vector_typename_size), sizeof(size_t));

  char read_vector_typename[vector_typename_size + 1];
  vecfp.read(reinterpret_cast<char *>(read_vector_typename), vector_typename_size);
  read_vector_typename[vector_typename_size] = '\0';

  if (std::strcmp(read_vector_typename, vector_typename.c_str()) != 0)
  {
    std::ostringstream err;
    err << "Wrong Vector format, found \'" << read_vector_typename << "\' while expecting \'"
        << vector_typename << "\'";
    print_error(err.str());
    exit(-1);
  }

  /* Read fileformat {BIN, SNAPPY} */
  int fileformat;
  vecfp.read(reinterpret_cast<char *>(&fileformat), sizeof(int));

  /* Vector dimension */
  size_t m;
  vecfp.read(reinterpret_cast<char *>(&m), sizeof(size_t));

  if (m * sizeof(ValueT) < props->ram_limit)
  {
    this->resize(m, 0);
  }
  else
  {
    print_error("Requested size is beyond \'ram_limit\'");
    vecfp.close();
    exit(-1);
  }

  if (fileformat == Utils::BIN)
  {
    vecfp.read(reinterpret_cast<char *>(this->data()), this->size() * sizeof(ValueT));
  }
  else if (fileformat == Utils::SNAPPY)
  {
    bool uncomp_succeed;

    size_t vec_snappy_size;
    vecfp.read(reinterpret_cast<char *>(&vec_snappy_size), sizeof(size_t));

    char *vec_snappy = new char[vec_snappy_size];
    vecfp.read(reinterpret_cast<char *>(vec_snappy), vec_snappy_size);

    uncomp_succeed = snappy::RawUncompress64(vec_snappy, vec_snappy_size, reinterpret_cast<char *>(this->data()));
    delete[] vec_snappy;

    if (!uncomp_succeed)
    {
      print_error("SNAPPY uncompression of vec failed");

      vecfp.close();
      exit(-1);
    }
  }

  vecfp.close();
}

template <typename ValueT>
Vector<ValueT> &Vector<ValueT>::operator=(Vector<ValueT> &&rvec)
{
  std::cout << "Move assignment" << std::endl;
  std::vector<ValueT>::operator = (std::move(rvec));
  std::swap(props, rvec.props);
  return *this;
}

template <typename ValueT>
Vector<ValueT> &Vector<ValueT>::operator+=(const Vector<ValueT>& rvec)
{
  if (this->size() != rvec.size())
  {
    std::ostringstream oss;
    oss << "Different sizes of left and right Vectors";
    print_error(oss.str());
    exit(-1);
  }

#pragma omp parallel for num_threads(props->nthreads)
  for (uint64_t i = 0; i < this->size(); i++)
  {
    this->at(i) += rvec[i];
  }

  return (*this);
}

template <typename ValueT>
Vector<ValueT> &Vector<ValueT>::operator+=(ValueT val)
{
#pragma omp parallel for num_threads(props->nthreads)
  for (uint64_t i = 0; i < this->size(); i++)
  {
    this->at(i) += val;
  }

  return (*this);
}

template <typename ValueT>
Vector<ValueT> &Vector<ValueT>::operator*=(ValueT val)
{
#pragma omp parallel for num_threads(props->nthreads)
  for (uint64_t i = 0; i < this->size(); i++)
  {
    this->at(i) *= val;
  }

  return (*this);
}

template <typename ValueT>
Vector<ValueT> &Vector<ValueT>::operator/=(Vector<ValueT> &rvec)
{
#pragma omp parallel for num_threads(props->nthreads)
  for (uint64_t i = 0; i < props->window; i++)
  {
    rvec[i] == 0 ? this->at(i) = 0 : this->at(i) /= rvec[i];
  }

  return (*this);
}

template <typename ValueT>
Vector<ValueT> &Vector<ValueT>::divide_and_sum_Nan(Vector<ValueT> &rvec, ValueT& aggs)
{
#pragma omp parallel for num_threads(props->nthreads) reduction(+:aggs)
  for (uint64_t i = 0; i < props->window; i++)
  {
    if(rvec[i]==0){
      aggs+=this->at(i);
      this->at(i) = 0;
    }else{
      this->at(i) /= rvec[i];
    }
  }

  return (*this);
}

template <typename ValueT>
Properties* Vector<ValueT>::get_properties()
{
  return props;
}

template <typename ValueT>
uint64_t Vector<ValueT>::get_lines() const
{
  return this->size();
}

} // namespace graphee

#endif // GRAPHEE_VECTOR_HPP__
