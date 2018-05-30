#ifndef GPE_VEC_H
#define GPE_VEC_H

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>

#include "gpe_props.h"
#include "gpe_utils.h"

namespace graphee {

template <typename val_t>
class gpe_vec : public std::vector<val_t> {
public:
  enum {BIN, SNAPPY};

  gpe_vec (gpe_props props) : 
    std::vector<val_t>(), props(props) {}

  gpe_vec (gpe_props props, size_t n, val_t init_val = 0) :
    std::vector<val_t> (n, init_val), props(props) {}

  ~gpe_vec () {}

  void save (std::string name, int fileformat = BIN, uint64_t offl = 0);
  void load (std::string name);

  template <typename gpe_mat_t>
  void mat_vec_prod (gpe_mat_t& mat, gpe_vec<val_t>& vec);

  typedef val_t value_type;

  const std::string vector_type {"GPE_VEC"};

private:
  std::ostringstream log;
  std::ostringstream wrn;
  std::ostringstream err;

  gpe_props props;
  uint64_t offl;

}; // class gpe_vec

template <class val_t>
void gpe_vec<val_t>::save (std::string name, int fileformat, uint64_t offl) {
  std::ofstream vecfp (name, std::ios_base::binary);

  size_t vector_type_size = vector_type.size();

  /* Save explicitly vector properties */
  vecfp.write (reinterpret_cast<const char*>(&vector_type_size), sizeof(size_t));
  vecfp.write (reinterpret_cast<const char*>(vector_type.c_str()), vector_type_size);

  /* Save fileformat {BIN, SNAPPY} */
  vecfp.write (reinterpret_cast<const char*>(&fileformat), sizeof(int));

  /* Printing offsets of the vector if needed */
  vecfp.write (reinterpret_cast<const char*>(&offl), sizeof(uint64_t));

  /* Vector dimension */
  size_t vec_size = this->size();
  vecfp.write (reinterpret_cast<const char*>(&vec_size), sizeof(size_t));

  if (fileformat == BIN) {
    vecfp.write (reinterpret_cast<const char*>(this->data()), this->size()*sizeof(val_t));

  } else if (fileformat == SNAPPY) {
    size_t vec_snappy_size = max_compress_size(this->size()*sizeof(val_t));
    char* vec_snappy = new char [vec_snappy_size];
    compress_snappy (reinterpret_cast<char*>(this->data()), this->size()*sizeof(val_t), vec_snappy, vec_snappy_size);

    vecfp.write (reinterpret_cast<const char*>(&vec_snappy_size), sizeof(size_t));
    vecfp.write (reinterpret_cast<const char*>(vec_snappy), vec_snappy_size);
    delete[] vec_snappy;
  }

  vecfp.close();
}

template <class val_t>
void gpe_vec<val_t>::load (std::string name) {
  std::ifstream vecfp (name, std::ios_base::binary);

  /* Save explicitly vector properties */
  size_t vector_type_size;
  vecfp.read (reinterpret_cast<char*>(&vector_type_size), sizeof(size_t));

  char read_vector_type[vector_type_size+1];
  vecfp.read (reinterpret_cast<char*>(read_vector_type), vector_type_size);
  read_vector_type[vector_type_size] = '\0';

  if (std::strcmp (read_vector_type, vector_type.c_str()) != 0) {
    err.str("");
    err << "Wrong vector format, found \'" << read_vector_type << "\' while expecting \'"
      << vector_type << "\'";
    gpe_error (err.str());
    exit (-1);
  }

  /* Read fileformat {BIN, SNAPPY} */
  int fileformat;
  vecfp.read (reinterpret_cast<char*>(&fileformat), sizeof(int));

  /* Reading offsets of the vector if needed */
  vecfp.read (reinterpret_cast<char*>(&offl), sizeof(uint64_t));

  /* Vector dimension */
  size_t n;
  vecfp.read (reinterpret_cast<char*>(&n), sizeof(size_t));

  if (n*sizeof(val_t) < props.ram_limit) {
    this->resize (n, 0);
  } else {
    gpe_error ("Requested size is beyond \'ram_limit\'");
    vecfp.close();
    exit (-1);
  }

  if (fileformat == BIN) {
    vecfp.read (reinterpret_cast<char*>(this->data()), this->size()*sizeof(val_t));

  } else if (fileformat == SNAPPY) {
    bool uncomp_succeed;

    size_t vec_snappy_size;
    vecfp.read (reinterpret_cast<char*>(&vec_snappy_size), sizeof(size_t));

    char* vec_snappy = new char [vec_snappy_size];
    vecfp.read (reinterpret_cast<char*>(vec_snappy), vec_snappy_size);

    uncomp_succeed = uncompress_snappy (vec_snappy, vec_snappy_size, reinterpret_cast<char*>(this->data()), n*sizeof(val_t));
    delete[] vec_snappy;

    if (!uncomp_succeed) {
      gpe_error("SNAPPY uncompression of vec failed");

      vecfp.close();
      exit (-1);
    }
  }
  
  vecfp.close();
}

template <typename val_t>
template <typename gpe_mat_t>
void gpe_vec<val_t>::mat_vec_prod (gpe_mat_t& mat, gpe_vec<val_t>& vec) {
  for (uint64_t id = 0; id < mat.m; id++) {
    for (uint64_t colid = mat.ia[id]; colid < mat.ia[id+1]; colid++) {
      this->at(id) += vec[mat.ja[colid]];
    }
  }
}

} // namespace graphee

#endif
