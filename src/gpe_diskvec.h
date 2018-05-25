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

template <class val_t>
class gpe_diskvec {
public:
  gpe_diskvec (gpe_props& i_props);
  gpe_diskvec (gpe_props& i_props, val_t init_val = 0);
  ~gpe_diskvec ();

private:
  std::fstream fp;

  std::vector <val_t> vec;

  void load_slice (uint64_t sliceID);
  void save_slice (uint64_t sliceID);

  void open_fp (uint64_t sliceID, std::ios_base::openmode mode);
  void close_fp (uint64_t sliceID);
}; // class gpe_diskvec

/*! Basic constructor of the class */
template <class val_t>
gpe_diskvec<val_t>::gpe_diskvec (gpe_props& i_props) {
  if (i_props.window*sizeof(val_t) > i_props.ram_limit) {
    gpe_error ("The \'gpe_vec\' size exceeds the \'ram_limit\'");
    exit (-1);
  }
  props = i_props;
}

/*! Complete constructor of the class */
template <class val_t>
gpe_diskvec<val_t>::gpe_diskvec (gpe_props& i_props, val_t init_val = 0) {
  if (i_props.window*sizeof(val_t) > i_props.ram_limit) {
    gpe_error ("The \'gpe_vec\' size exceeds the \'ram_limit\'");
    exit (-1);
  }
  props = i_props;
  /*
   * HERE COMPLETE THE FULLFILL OF THE VECTOR !
   */
}

} // namesapce graphee

#endif // GPE_DISKVEC_H
