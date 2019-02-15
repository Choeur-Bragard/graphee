#ifndef GPE_PROPERTIES_HPP__
#define GPE_PROPERTIES_HPP__

#include <atomic>
#include <iostream>
#include <string>

namespace graphee {

class Properties {
public:
  Properties()
      : name(""), nvertices(0), declared_nvertices(0), nslices(0), nthreads(0),
        ram_limit(0), sort_limit(0), nblocks(0), window(0) {}

  Properties(std::string name, uint64_t declared_nvertices, uint64_t nslices,
             uint64_t nthreads, size_t ram_limit, size_t sort_limit)
      : name(name), nvertices(declared_nvertices % nslices == 0
                                  ? declared_nvertices
                                  : declared_nvertices + nslices -
                                        declared_nvertices % nslices),
        declared_nvertices(declared_nvertices), nslices(nslices),
        nthreads(nthreads), ram_limit(ram_limit), sort_limit(sort_limit),
        nblocks(nslices * nslices), window(nvertices / nslices),
        alloc_memory(0) {}

  ~Properties() {}

  static const size_t B{1};
  static const size_t KB{1UL << 10};
  static const size_t MB{1UL << 20};
  static const size_t GB{1UL << 30};

  const std::string name;
  // the matrix size. (can be bigger than the declared number of
  // vertice in order to have a multiple of nslices)
  const uint64_t nvertices;
  // this is the actual number of vertices in the graph
  const uint64_t declared_nvertices;
  const uint64_t nslices;
  const uint64_t nthreads;
  const size_t ram_limit;
  const size_t sort_limit;
  const uint64_t nblocks;
  const uint64_t window;

  std::atomic<size_t> alloc_memory;
}; // class graphee::Properties

} // namespace graphee

#endif // GPE_PROPERTIES_HPP__
