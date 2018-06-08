#ifndef GRAPHEE_VECTOR_H__
#define GRAPHEE_VECTOR_H__

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>

#include "utils.h"
#include "properties.h"

#include "sparse_matrix_csr.h"

namespace graphee
{

/*! \brief Graphee vector class
 *         Includes load/dump of data on disk
 *
 * It contains overloaded operators for operations
 * with `sparseMatrix*` classes
 */

template <typename valueT>
class vector : public std::vector<valueT>
{
public:
  vector() : std::vector<valueT>() {}
  vector(properties &properties, size_t m, valueT init_value = 0) : std::vector<valueT>(m, init_value), m(m), props(properties) {}
  ~vector() {}

  void save(std::string name, int fileformat = utils::BIN);
  void load(std::string name);

  vector<valueT> &operator+=(vector<valueT> rvec);
  vector<valueT> &operator+=(valueT val);

  const std::string vectorType{"vector"};

  using valueType = valueT;

private:
  properties &&props;
  uint64_t m;

}; // class graphee::vector

template <class valueT>
void vector<valueT>::save(std::string name, int fileformat)
{
  std::ofstream vecfp(name, std::ios_base::binary);

  size_t vectorType_size = vectorType.size();

  /* Save explicitly vector properties */
  vecfp.write(reinterpret_cast<const char *>(&vectorType_size), sizeof(size_t));
  vecfp.write(reinterpret_cast<const char *>(vectorType.c_str()), vectorType_size);

  /* Save fileformat {BIN, SNAPPY} */
  vecfp.write(reinterpret_cast<const char *>(&fileformat), sizeof(int));

  /* Vector dimension */
  size_t vec_size = this->size();
  vecfp.write(reinterpret_cast<const char *>(&vec_size), sizeof(size_t));

  if (fileformat == utils::BIN)
  {
    vecfp.write(reinterpret_cast<const char *>(this->data()), this->size() * sizeof(valueT));
  }
  else if (fileformat == utils::SNAPPY)
  {
    size_t vec_snappy_size = max_compress_size(this->size() * sizeof(valueT));
    char *vec_snappy = new char[vec_snappy_size];
    snappy::RawCompress64(reinterpret_cast<char *>(this->data()), this->size() * sizeof(valueT), vec_snappy, vec_snappy_size);

    vecfp.write(reinterpret_cast<const char *>(&vec_snappy_size), sizeof(size_t));
    vecfp.write(reinterpret_cast<const char *>(vec_snappy), vec_snappy_size);
    delete[] vec_snappy;
  }

  vecfp.close();
}

template <class valueT>
void vector<valueT>::load(std::string name)
{
  this->clear();
  std::ifstream vecfp(name, std::ios_base::binary);

  /* Save explicitly vector properties */
  size_t vectorType_size;
  vecfp.read(reinterpret_cast<char *>(&vectorType_size), sizeof(size_t));

  char read_vectorType[vectorType_size + 1];
  vecfp.read(reinterpret_cast<char *>(read_vectorType), vectorType_size);
  read_vectorType[vectorType_size] = '\0';

  if (std::strcmp(read_vectorType, vectorType.c_str()) != 0)
  {
    std::ostringstream err;
    err << "Wrong vector format, found \'" << read_vectorType << "\' while expecting \'"
        << vectorType << "\'";
    print_error(err.str());
    exit(-1);
  }

  /* Read fileformat {BIN, SNAPPY} */
  int fileformat;
  vecfp.read(reinterpret_cast<char *>(&fileformat), sizeof(int));

  /* Vector dimension */
  size_t m;
  vecfp.read(reinterpret_cast<char *>(&m), sizeof(size_t));

  if (m * sizeof(valueT) < props.ram_limit)
  {
    this->resize(m, 0);
  }
  else
  {
    print_error("Requested size is beyond \'ram_limit\'");
    vecfp.close();
    exit(-1);
  }

  if (fileformat == utils::BIN)
  {
    vecfp.read(reinterpret_cast<char *>(this->data()), this->size() * sizeof(valueT));
  }
  else if (fileformat == utils::SNAPPY)
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

template <typename valueT>
vector<valueT> &vector<valueT>::operator+=(vector<valueT> rvec)
{
  if (this->size() != rvec.size())
  {
    std::ostringstream oss;
    oss << "Different sizes of left and right vectors";
    print_error(oss.str());
    exit(-1);
  }

  for (uint64_t i = 0; i < this->size(); i++)
  {
    this->at(i) += rvec[i];
  }

  return (*this);
}

template <typename valueT>
vector<valueT> &vector<valueT>::operator+=(valueT val)
{
  for (uint64_t i = 0; i < this->size(); i++)
  {
    this->at(i) += val;
  }

  return (*this);
}

} // namespace graphee

#endif // GRAPHEE_VECTOR_H__
