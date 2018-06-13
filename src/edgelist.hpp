#ifndef GRAPHEE_EDGELIST_HPP__
#define GRAPHEE_EDGELIST_HPP__

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <mutex>
#include <thread>
#include <array>

#include <cstdio>
#include <cassert>
#include <zlib.h>

#include "properties.hpp"
#include "utils.hpp"

namespace graphee
{

class Edgelist
{
public:
  Edgelist (Properties *properties, std::mutex *mutex,
            size_t sstream_limit, std::vector<std::string> &filelist) :
    props(properties), mtx(mutex), sstream_limit(sstream_limit),
    filelist(std::move(filelist)), file_ptr(nullptr) {}

  bool read_file(std::stringstream& sstream_read);

private:
  Properties *props;
  std::mutex *mtx;

  std::vector<std::string> filelist;
  std::vector<std::string>::iterator filelist_it;

  const size_t sstream_limit;

  FILE* file_ptr;

  static void deflate_chunk(Edgelist* edgelist, std::stringstream &sstream_read);
}; // class Edgelist

} // namespace graphee

#endif // GRAPHEE_EDGELIST_HPP__
