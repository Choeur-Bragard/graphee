#ifndef GRAPHEE_DISKVEC_H__
#define GRAPHEE_DISKVEC_H__

#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

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

namespace graphee
{

template <typename vectorT>
class disk_vector
{
public:
  disk_vector() {}
  disk_vector(properties &properties,
              std::string vector_name,
              typename vectorT::valueType init_val = 0)
      : props(properties), name(vector_name)
  {
    if (props.window * sizeof(typename vectorT::valueType) > props.ram_limit)
    {
      print_error("The \'graphee::vector\' size exceeds the \'ram_limit\'");
      exit(-1);
    }

    vectorT tmp(props, props.window, init_val);
    for (uint64_t slice_id = 0; slice_id < props.nslices; i++)
    {
      tmp.save(get_slice_filename(slice_id));
    }
  }

  ~disk_vector() {}

  vectorT &get_slice(uint64_t slice_id);

  void swap(disk_vector<vectorT> &rvec);

  disk_vector<vectorT> &operator+=(typename vectorT::valueType val);

  using vectorType = vectorT;

private:
  properties &&props;
  std::string name;

  std::string get_slice_filename(uint64_t slice_id);
}; // class disk_vector

template <typename vectorT>
vectorT &
disk_vector<vectorT>::get_slice(uint64_t slice_id)
{
  std::ostringstream oss;
  oss << "Load slice [" << slice_id << "] of vector \'" << name << "\'";
  print_log(oss.str());

  vectorT res;
  res.load(get_slice_filename(slice_id));
}

template <typename vectorT>
std::string
disk_vector<vectorT>::get_slice_filename(uint64_t slice_id)
{
  std::ostringstream slicename;
  slicename << props.name << "_" << vec_name << "_dvecslc_" << slice_id
            << ".gpe";
  return slicename.str();
}

template <typename vectorT>
void disk_vector<vectorT>::swap(disk_vector<vectorT> &vec)
{
  if (props.nvertices != vec.props.nvertices)
  {
    err.str("");
    err << "Could not swap vector \'" << vec_name << "\' and \'";
    err << vec.vec_name << "\' because dimensions are not equal";
    gpe_error(err.str());
  }

  int succed;
  std::ostringstream tmpname;
  tmpname << vec_name << "_swap_file.gpe";
  for (uint64_t sliceID = 0; sliceID < props.nslices; sliceID++)
  {
    succed = rename(get_slice_filename(sliceID).c_str(), tmpname.str().c_str());
    if (succed != 0)
    {
      err.str("");
      err << "Could not swap vector \'" << get_slice_filename(sliceID)
          << "\' to \'";
      err << tmpname.str() << "\'";
      gpe_error(err.str());
    }

    succed = rename(vec.get_slice_filename(sliceID).c_str(),
                    get_slice_filename(sliceID).c_str());
    if (succed != 0)
    {
      err.str("");
      err << "Could not swap vector \'" << vec.get_slice_filename(sliceID)
          << "\' to \'";
      err << get_slice_filename(sliceID) << "\'";
      gpe_error(err.str());
    }

    succed =
        rename(tmpname.str().c_str(), vec.get_slice_filename(sliceID).c_str());
    if (succed != 0)
    {
      err.str("");
      err << "Could not swap vector \'" << tmpname.str() << "\' to \'";
      err << vec.get_slice_filename(sliceID) << "\'";
      gpe_error(err.str());
    }
  }
}

template <typename vectorT>
template <typename gpe_dmat_t, typename gpe_dvec_t>
void disk_vector<vectorT>::alpha_mat_vec_prod(float alpha,
                                              gpe_dmat_t &dmat,
                                              gpe_dvec_t &dvec)
{
  vectorT res(props);
  typename gpe_dvec_t::vector_type vec_arg(props);
  typename gpe_dmat_t::matrix_type mat_arg(props);

  for (uint64_t line = 0; line < props.nslices; line++)
  {
    get_vector_slice(line, res);
    for (uint64_t col = 0; col < props.nslices; col++)
    {
      dvec.get_vector_slice(col, vec_arg);
      dmat.get_matrix_block(line, col, mat_arg);

      res.alpha_mat_vec_prod(alpha, mat_arg, vec_arg);
    }
    res.save(get_slice_filename(line));
  }
}

template <typename vectorT>
void disk_vector<vectorT>::operator+=(typename vectorT::value_type val)
{
  vectorT vec(props);
  for (uint64_t slice_id = 0; slice_id < props.nslices; slice_id++)
  {
    vec.load(get_slice_filename(slice_id));
    vec += val;
    vec.save(get_slice_filename(slice_id));
  }
}

template <typename vectorT>
void disk_vector<vectorT>::set_init_value(typename vectorT::value_type init_val)
{
  log.str("");
  log << "Init value of vector \'" << vec_name << "\'";
  gpe_log(log.str());
  for (uint64_t sliceID = 0; sliceID < props.nslices; sliceID++)
  {
    vectorT vec(props, props.window, init_val);
    vec.save(get_slice_filename(sliceID));
  }
}

} // namespace graphee

#endif // GRAPHEE_DISKVEC_H__
