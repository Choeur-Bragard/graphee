#ifndef GPE_DISKVEC_H
#define GPE_DISKVEC_H

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>

#include "gpe_utils.h"
#include "gpe_diskmat.h"
#include "gpe_vec.h"

/* \brief Vector saved by slices within disk
 *
 * Modify the format in order to template a vector 
 * standart class instead of a value type.
 * In such cases either sparse or dense vector could be written,
 * on the disk.
 *
 */

namespace graphee {

template <class gpe_vec_t>
class gpe_diskvec {
public:
  gpe_diskvec (gpe_props props, std::string vec_name); 
  gpe_diskvec (gpe_props props, std::string vec_name, size_t n, 
    typename gpe_vec_t::value_type init_val = 0);

  ~gpe_diskvec () {}

  void get_vector_slice (uint64_t sliceID, gpe_vec_t& vec);

private:
  gpe_props props;
  std::string vec_name;

  std::ostringstream log;
  std::ostringstream wrn;
  std::ostringstream err;

  std::string get_slice_filename (uint64_t sliceID);
}; // class gpe_diskvec

template <class gpe_vec_t>
gpe_diskvec<gpe_vec_t>::gpe_diskvec (gpe_props arg_props, std::string arg_vec_name) {
  props = arg_props;
  vec_name = arg_vec_name;

  if (props.window*sizeof(typename gpe_vec_t::value_type) > props.ram_limit) {
    gpe_error ("The \'gpe_vec\' size exceeds the \'ram_limit\'");
    exit (-1);
  }
}

template <class gpe_vec_t>
gpe_diskvec<gpe_vec_t>::gpe_diskvec (gpe_props arg_props, std::string arg_vec_name, size_t n, 
  typename gpe_vec_t::value_type init_val) {
  props = arg_props;
  vec_name = arg_vec_name;

  if (props.window*sizeof(typename gpe_vec_t::value_type) > props.ram_limit) {
    gpe_error ("The \'gpe_vec\' size exceeds the \'ram_limit\'");
    exit (-1);
  }

  for (uint64_t sliceID = 0; sliceID < props.nslices; sliceID++) {
    gpe_vec_t vec (props, props.window, init_val);
    vec.save (get_slice_filename(sliceID));
  }
}

template <class gpe_vec_t>
void gpe_diskvec<gpe_vec_t>::get_vector_slice (uint64_t sliceID, gpe_vec_t& vec) {
  log.str("");
  log << "Start to load vector slice [" << sliceID << "]";
  gpe_log (log.str());

  vec.load (get_slice_filename(sliceID));

  log.str("");
  log << "Loaded vector slice [" << sliceID << "]";
  gpe_log (log.str());
}

template <class gpe_vec_t>
std::string gpe_diskvec<gpe_vec_t>::get_slice_filename (uint64_t sliceID) {
  std::ostringstream slicename;
  slicename << props.name << "_" << vec_name << "_dvecslc_" << sliceID << ".gpe";
  return slicename.str();
}

} // namesapce graphee

#endif // GPE_DISKVEC_H
