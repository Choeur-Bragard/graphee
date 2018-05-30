#ifndef GPE_DMAT_DVEC_MULT_H
#define GPE_DMAT_DVEC_MULT_H

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>

#include "graphee.h"

namespace graphee {

template <typename gpe_dmat_t, typename gpe_dvec_t>
class gpe_dmat_dvec_mult {
public:
  gpe_dmat_dvec_mult (gpe_props props);
  ~gpe_dmat_dvec_mult () {}
}; // class gpe_dmat_dvec_mult

} // namespace graphee

#endif //GPE_DMAT_DVEC_MULT_H
