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
  disk_vector(properties &properties, std::string vector_name,
              typename vectorT::valueType init_val = 0)
    : props(properties), name(vector_name)
  {
    if (props.window * sizeof(typename vectorT::valueType) > props.ram_limit)
    {
      print_error("The \'graphee::vector\' size exceeds the \'ram_limit\'");
      exit(-1);
    }

    vectorT tmp(props, props.window, init_val);
    for (uint64_t slice_id = 0; slice_id < props.nslices; slice_id++)
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
};  // class disk_vector

template <typename vectorT>
vectorT &disk_vector<vectorT>::get_slice(uint64_t slice_id)
{
  std::ostringstream oss;
  oss << "Load slice [" << slice_id << "] of vector \'" << name << "\'";
  print_log(oss.str());

  vectorT res;
  res.load(get_slice_filename(slice_id));
}

template <typename vectorT>
std::string disk_vector<vectorT>::get_slice_filename(uint64_t slice_id)
{
  std::ostringstream slicename;
  slicename << props.name << "_" << name << "_dvecslc_" << slice_id
            << ".gpe";
  return slicename.str();
}

template <typename vectorT>
void disk_vector<vectorT>::swap(disk_vector<vectorT> &vec)
{
  if (props.nvertices != vec.props.nvertices)
  {
    std::ostringstream oss;
    oss << "Could not swap vector \'" << name << "\' and \'";
    oss << vec.name << "\' because dimensions are not equal";
    print_error(oss.str());
  }

  int succed;
  std::ostringstream tmpname;
  tmpname << name << "_swap_file.gpe";
  for (uint64_t slice_id = 0; slice_id < props.nslices; slice_id++)
  {
    succed = rename(get_slice_filename(slice_id).c_str(),
                    tmpname.str().c_str());

    if (succed != 0)
    {
      std::ostringstream oss;
      oss << "Could not swap vector \'" << get_slice_filename(slice_id)
          << "\' to \'";
      oss << tmpname.str() << "\'";
      print_error(oss.str());
    }

    succed = rename(vec.get_slice_filename(slice_id).c_str(),
                    get_slice_filename(slice_id).c_str());

    if (succed != 0)
    {
      std::ostringstream oss;
      oss << "Could not swap vector \'" << vec.get_slice_filename(slice_id)
          << "\' to \'";
      oss << get_slice_filename(slice_id) << "\'";
      print_error(oss.str());
    }

    succed = rename(tmpname.str().c_str(),
                    vec.get_slice_filename(slice_id).c_str());

    if (succed != 0)
    {
      std::ostringstream oss;
      oss << "Could not swap vector \'" << tmpname.str() << "\' to \'";
      oss << vec.get_slice_filename(slice_id) << "\'";
      print_error(oss.str());
    }
  }
}

template <typename vectorT>
disk_vector<vectorT>& disk_vector<vectorT>::operator+=(typename vectorT::valueType val)
{
  vectorT vec(props);
  for (uint64_t slice_id = 0; slice_id < props.nslices; slice_id++)
  {
    vec.load(get_slice_filename(slice_id));
    vec += val;
    vec.save(get_slice_filename(slice_id));
  }

  return (*this);
}

}  // namespace graphee

#endif  // GRAPHEE_DISKVEC_H__
