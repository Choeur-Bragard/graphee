
#include "graphee.hpp"

#include <iostream>
#include <string>

/** Dump a pagerank vector
 * usage ./dump_pr vector_part0.gpe vector_part1.gpe ...
 */
int main(int argc, char **argv) {

  graphee::Properties props(
      std::string("dump"), // name of your graph
      0,                           // number of nodes
      5,                                // number of slices
      4,                                // number of threads
      5 * graphee::Properties::GB,      // max RAM value
      32 * graphee::Properties::MB);    // max size of sorting vector

  graphee::Vector<double> vec(&props);

  uint64_t n = 0;
  double sum = 0;

  for (int i = 1; i < argc; i++) {
    std::string file(argv[i]);
    vec.load(file);

    for (double v : vec) {
      std::cout << n << "\t" << v << std::endl;
      n++;
      sum += v;
    }
  }
}