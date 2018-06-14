#include "edgelist.hpp"

namespace graphee
{

bool Edgelist::read(std::stringstream &sstream_read)
{
  mtx->lock();

  sstream_read.str("");
  sstream_read.clear();

  bool has_file;

  // C style code zone, NEEDED FOR ZLIB
  if (file_ptr == nullptr) // not initialized yet
  {
    filelist_it = filelist.begin();
    file_ptr = gzopen(filelist_it->c_str(), "rb");
    part_id = 0;
    has_file = true;
  }

  if (gzeof(file_ptr))
  {
    gzclose(file_ptr);
    filelist_it++;
    if (filelist_it != filelist.end())
    {
      file_ptr = gzopen(filelist_it->c_str(), "rb");
      part_id = 0;
      has_file = true;
    }
    else
    {
      has_file = false;
    }
  }
  
  if (file_ptr == Z_NULL)
  {
    std::ostringstream oss;
    oss << "Cannot open file \'" << *filelist_it << "\'";
    print_error(oss.str());
    exit(-1);
  }

  // End of C style code zone

  if (has_file)
  {
    part_id++;
    std::thread t (deflate_chunk, this, std::ref(sstream_read));
    t.detach();
  }

  return has_file || buf_filed; 
}

std::string Edgelist::get_filename()
{
  return *filelist_it;
}

void Edgelist::deflate_chunk(Edgelist* el, std::stringstream &sstream_read)
{
  std::ostringstream oss;
  oss << "Defalte part " << el->part_id << " of file \'" << *(el->filelist_it) << "\'";
  print_log(oss.str());

  int ret = gzread(el->file_ptr, el->buf, el->buf_size);

  if (ret > 0)
    sstream_read << el->buf;

  el->buf_filed = ret > 0;

  el->mtx->unlock();
}

} // namespace graphee
