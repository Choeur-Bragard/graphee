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

template <class val_t>
class gpe_vec : public std::vector<val_t> {
public:
  enum {BIN, SNAPPY};

  gpe_vec (gpe_props& i_prop) : 
    std::vector<val_t>() {}

  gpe_vec (gpe_props& i_prop, size_t n, val_t init_val = 0) :
    std::vector<val_t> (n, init_val) {}

  ~gpe_vec ();

  void save (std::string name, int fileformat = BIN, uint64_t offl = 0);
  void load (std::string name);

  const typedef val_t value_type;
  const std::string vector_type {"GPE_VEC"};

private:
  uint64_t offl;

}; // class gpe_vec

template <class val_t>
void gpe_vec<val_t>::save (std::string name, int fileformat = BIN, uint64_t offl = 0) {
  std::ofstream vecfp (name, std::ios_base::binary);

  size_t vec_type_size = vec_type.size();

  /* Save explicitly vector properties */
  vecfp.write (reinterpret_cast<const char*>(&vec_type_size), sizeof(size_t));
  vecfp.write (reinterpret_cast<const char*>(vector_type.c_str()), vector_type_size);

  /* Save fileformat {BIN, SNAPPY} */
  vecfp.write (reinterpret_cast<const char*>(&fileformat), sizeof(int));

  /* Printing offsets of the vector if needed */
  vecfp.write (reinterpret_cast<const char*>(&offl), sizeof(uint64_t));

  /* Vector dimension */
  vecfp.write (reinterpret_cast<const char*>(&(this->size())), sizeof(size_t));

  if (fileformat == BIN) {
    vecfp.write (reinterpret_cast<const char*>(this->data()), this->size()*sizeof(val_t));

  } else if (fileformat == SNAPPY) {
    char* vec_snappy = new char [snappy::MaxCompressedLength(this->size()*sizeof(val_t))];
    size_t vec_snappy_size;
    compress_snappy (reinterpret_cast<char*>(this->data()), this->size()*sizeof(val_t), vec_snappy, vec_snappy_size);

    vecfp.write (reinterpret_cast<const char*>(&vec_snappy_size), sizeof(size_t));
    vecfp.write (reinterpret_cast<const char*>(vec_snappy), vec_snappy_size);
    delete[] vec_snappy;
  }

  vec.close();
}

template <class val_t>
void gpe_diskvec<val_t>::load (std::string name) {
  std::ifstream vecfp (name, std::ios_base::binary);

  /* Save explicitly vector properties */
  size_t vec_type_size;
  vecfp.read (reinterpret_cast<char*>(&vec_type_size), sizeof(size_t));

  char read_vec_type[vec_type_size];
  vecfp.read (reinterpret_cast<char*>(read_vector_type), vec_type_size);

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
    matfp.read (reinterpret_cast<char*>(this->data()), this->size()*sizeof(val_t));

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

} // namespace graphee

#endif
