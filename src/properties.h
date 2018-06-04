#ifndef GPE_PROPS_H
#define GPE_PROPS_H

#include <iostream>
#include <string>

namespace graphee {

class properties {
public:
  properties () : name(""), nvertices(0), nslices(0), nthreads(0), ram_limit(0), sort_limit(0),
    nblocks(0), window(0) {}

  properties (std::string name, uint64_t nvertices, uint64_t nslices, uint64_t nthreads, 
      size_t ram_limit, size_t sort_limit) :

    name(name), nvertices(nvertices), nslices(nslices), nthreads(nthreads),
    ram_limit(ram_limit), sort_limit(sort_limit), nblocks(nslices*nslices), 
    window(nvertices/nslices) {}

  ~gpe_props () {}

  const size_t B  {1};
  const size_t KB {1UL << 10};
  const size_t MB {1UL << 20};
  const size_t GB {1UL << 30};

  const std::string name;
  const uint64_t nvertices;
  const uint64_t nslices;
  const uint64_t nthreads;
  const size_t ram_limit;
  const size_t sort_limit;
  const uint64_t nblocks;
  const uint64_t window;
};

} // namespace graphee

#endif // SETINGS_H
