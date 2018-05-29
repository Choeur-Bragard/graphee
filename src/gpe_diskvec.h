#ifndef GPE_DISKVEC_H
#define GPE_DISKVEC_H

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>

#include "gpe_utils.h"
#include "gpe_diskmat.h"

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
  gpe_diskvec (gpe_props props, std::string vec_name) :
    props(props), vec_name(vec_name) {}

  gpe_diskvec (gpe_props props, std::string vec_name size_t n, 
      gpe_vec_t::value_type init_val = 0) : gpe_diskvec (i_props, vec_name) {}

  ~gpe_diskvec ();

  void get_vector_block (uint64_t sliceID, gpe_vec_t& vec);

private:
  gpe_props props;
  std::string vec_name;

  std::string get_slice_filename (uint64_t sliceID);
}; // class gpe_diskvec

/*! Basic constructor of the class */
template <class gpe_vec_t>
gpe_diskvec<val_t>::gpe_diskvec (gpe_props props) {
  if (i_props.window*sizeof(gpe_vec_t::value_type) > i_props.ram_limit) {
    gpe_error ("The \'gpe_vec\' size exceeds the \'ram_limit\'");
    exit (-1);
  }
}

/*! Complete constructor of the class */
template <class gpe_vec_t>
gpe_diskvec<gpe_vec_t>::gpe_diskvec (gpe_props props, val_t init_val = 0) {
  for (uint64_t sliceID = 0; sliceID < props.nslices; sliceID++) {
    gpe_vec_t vec (props, props.window, init_val);
    vec.save (get_slice_filename(sliceID));
  }
}

} // namesapce graphee

#endif // GPE_DISKVEC_H
