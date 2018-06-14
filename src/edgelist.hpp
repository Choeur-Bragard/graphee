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
            const size_t buf_size, std::vector<std::string> &filelist) :
    props(properties), mtx(mutex), buf_size(buf_size), part_id(0),
    buf_filed(true), filelist(std::move(filelist)), file_ptr(nullptr)
  {
    buf = new char[buf_size];
  }

  bool read(std::stringstream& sstream_read);

  std::string get_filename();

private:
  Properties *props;
  std::mutex *mtx;

  std::vector<std::string> filelist;
  std::vector<std::string>::iterator filelist_it;

  size_t part_id;
  char* buf;
  const size_t buf_size;
  bool buf_filed;

  gzFile file_ptr;

  static void deflate_chunk(Edgelist* edgelist, std::stringstream &sstream_read);
}; // class Edgelist

} // namespace graphee

#endif // GRAPHEE_EDGELIST_HPP__
