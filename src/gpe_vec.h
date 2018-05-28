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

  gpe_vec (gpe_props& i_prop, std::string i_name) : 
    std::vector<val_t>() {}

  gpe_vec (gpe_props& i_prop, std::string i_name, size_t n, val_t init_val = 0) :
    std::vector<val_t> (n, init_val) {}

  ~gpe_vec ();

  void save (int fileformat = BIN);
  void load ();

private:
  std::string name;
}; // class gpe_vec

} // namespace graphee

#endif
