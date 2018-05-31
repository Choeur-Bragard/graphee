#ifndef GPE_PROPS_H
#define GPE_PROPS_H

#include <iostream>
#include <string>

namespace graphee {

class gpe_props {
public:
  enum {B, KB, MB, GB};

  gpe_props () {}

  gpe_props (std::string in_name, uint64_t in_nvertices, uint64_t in_nslices, uint64_t in_nthreads, 
      uint64_t in_ram_limit, int unitRL, uint64_t in_sort_limit, int unitSL) :

    name(in_name), nvertices(in_nvertices), nslices(in_nslices),
    nthreads(in_nthreads), ram_limit(in_ram_limit), sort_limit(in_sort_limit) {

      if (unitRL == KB) {
        ram_limit *= 1 << 10;
      } else if (unitRL == MB) {
        ram_limit *= 1 << 20;
      } else if (unitRL == GB) {
        ram_limit *= 1 << 30;
      }

      if (unitSL == KB) {
        sort_limit *= 1 << 10;
      } else if (unitSL == MB) {
        sort_limit *= 1 << 20;
      } else if (unitSL == GB) {
        sort_limit *= 1 << 30;
      }

      nvertices *= 1000000;
      window = nvertices/nslices;
      nblocks = nslices*nslices;
  }

  ~gpe_props () {}

  std::string name;
  uint64_t nthreads {0};
  uint64_t nvertices {0};
  uint64_t ram_limit {0};
  uint64_t sort_limit {0};
  uint64_t nslices {0};
  uint64_t nblocks {0};
  uint64_t window {0};
};

} // namespace graphee

#endif // SETINGS_H
