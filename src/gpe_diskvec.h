#ifndef GPE_DISKVEC_H
#define GPE_DISKVEC_H

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>

#include <cstdio>

#include "graphee.h"

/* \brief Vector saved by slices within disk
 *
 * Modify the format in order to template a vector 
 * standart class instead of a value type.
 * In such cases either sparse or dense vector could be written,
 * on the disk.
 *
 */

namespace graphee {

template <typename gpe_vec_t>
class gpe_diskvec {
public:
  gpe_diskvec (gpe_props props, std::string vec_name); 
  gpe_diskvec (gpe_props props, std::string vec_name, size_t n, 
    typename gpe_vec_t::value_type init_val = 0);

  ~gpe_diskvec () {}

  void get_vector_slice (uint64_t sliceID, gpe_vec_t& vec);

  void swap (gpe_diskvec<gpe_vec_t>& vec);

  template <typename gpe_dmat_t, typename gpe_dvec_t>
  void mat_vec_prod (gpe_dmat_t& dmat, gpe_dvec_t& dvec);

  typedef gpe_vec_t vector_type;

private:
  gpe_props props;
  std::string vec_name;

  std::ostringstream log;
  std::ostringstream wrn;
  std::ostringstream err;

  std::string get_slice_filename (uint64_t sliceID);
}; // class gpe_diskvec

template <typename gpe_vec_t>
gpe_diskvec<gpe_vec_t>::gpe_diskvec (gpe_props arg_props, std::string arg_vec_name) {
  props = arg_props;
  vec_name = arg_vec_name;

  if (props.window*sizeof(typename gpe_vec_t::value_type) > props.ram_limit) {
    gpe_error ("The \'gpe_vec\' size exceeds the \'ram_limit\'");
    exit (-1);
  }
}

template <typename gpe_vec_t>
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

template <typename gpe_vec_t>
void gpe_diskvec<gpe_vec_t>::get_vector_slice (uint64_t sliceID, gpe_vec_t& vec) {
  log.str("");
  log << "Start to load vector slice [" << sliceID << "]";
  gpe_log (log.str());

  vec.clear ();
  vec.load (get_slice_filename(sliceID));

  log.str("");
  log << "Loaded vector slice [" << sliceID << "]";
  gpe_log (log.str());
}

template <typename gpe_vec_t>
std::string gpe_diskvec<gpe_vec_t>::get_slice_filename (uint64_t sliceID) {
  std::ostringstream slicename;
  slicename << props.name << "_" << vec_name << "_dvecslc_" << sliceID << ".gpe";
  return slicename.str();
}

template <typename gpe_vec_t>
void gpe_diskvec<gpe_vec_t>::swap (gpe_diskvec<gpe_vec_t>& vec) {
  if (props.nvertices != vec.props.nvertices) {
    err.str("");
    err << "Could not swap vector \'" << vec_name << "\' and \'";
    err << vec.vec_name << "\' because dimensions are not equal";
    gpe_error (err.str());
  }

  int succed;
  std::ostringstream tmpname;
  tmpname << vec_name << "_swap_file.gpe";
  for (uint64_t sliceID = 0; sliceID < props.nslices; sliceID++) {
    succed = rename (get_slice_filename(sliceID).c_str(), tmpname.str().c_str());
    if (succed != 0) {
      err.str("");
      err << "Could not swap vector \'" << get_slice_filename(sliceID) << "\' to \'";
      err << tmpname.str() << "\'";
      gpe_error (err.str());
    }

    succed = rename (vec.get_slice_filename(sliceID).c_str(), get_slice_filename(sliceID).c_str());
    if (succed != 0) {
      err.str("");
      err << "Could not swap vector \'" << vec.get_slice_filename(sliceID) << "\' to \'";
      err << get_slice_filename(sliceID) << "\'";
      gpe_error (err.str());
    }

    succed = rename (tmpname.str().c_str(), vec.get_slice_filename(sliceID).c_str());
    if (succed != 0) {
      err.str("");
      err << "Could not swap vector \'" << tmpname.str() << "\' to \'";
      err << vec.get_slice_filename(sliceID) << "\'";
      gpe_error (err.str());
    }
  }
}

template <typename gpe_vec_t>
template <typename gpe_dmat_t, typename gpe_dvec_t>
void gpe_diskvec<gpe_vec_t>::mat_vec_prod (gpe_dmat_t& dmat, gpe_dvec_t& dvec) {
  gpe_vec_t res (props);
  typename gpe_dvec_t::vector_type vec_arg (props);
  typename gpe_dmat_t::matrix_type mat_arg (props);

  for (uint64_t line = 0; line < props.nslices; line++) {
    get_vector_slice (line, res);
    for (uint64_t col = 0; col < props.nslices; col++) {
      dvec.get_vector_slice (col, vec_arg);
      dmat.get_matrix_block (line, col, mat_arg);

      res.mat_vec_prod(mat_arg, vec_arg);
    }
    res.save (get_slice_filename(line));
  }
}

} // namesapce graphee

#endif // GPE_DISKVEC_H
