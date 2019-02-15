#include "edgelist.hpp"

namespace graphee
{

bool Edgelist::read(std::stringstream &sstream_read)
{
  mtx->lock();

  sstream_read.str("");
  sstream_read.clear();

  // C style code zone, NEEDED FOR ZLIB
  if (file_ptr == nullptr) // not initialized yet
  {
    filelist_it = filelist.begin();
    file_ptr = gzopen(filelist_it->c_str(), "rb");
    part_id = 0;
  }

  if (filelist_it != filelist.end())
  {
    if (gzeof(file_ptr))
    {
      gzclose(file_ptr);
      filelist_it++;
      if (filelist_it != filelist.end())
      {
        file_ptr = gzopen(filelist_it->c_str(), "rb");
        part_id = 0;
      }
      else
      {
        mtx->unlock();
        return is_data;
      }
    }

    if (file_ptr == Z_NULL)
    {
      std::ostringstream oss;
      oss << "Cannot open file \'" << *filelist_it << "\'";
      print_error(oss.str());
      exit(-1);
    }
    else
    {
      part_id++;
      std::thread t (deflate_chunk, this, std::ref(sstream_read));
      t.detach();
    }

    return true;
  }
  else
  {
    mtx->unlock();
    return false;
  }
}

std::string Edgelist::get_filename()
{
  return (filelist_it != filelist.end()) ? *filelist_it : *(filelist_it-1);
}

void Edgelist::deflate_chunk(Edgelist* el, std::stringstream &sstream_read)
{
  int ret = gzread(el->file_ptr, el->buf, el->buf_size);

  if (ret > 0)
    sstream_read << el->buf;

  std::ostringstream oss;
  oss << "Deflate part " << el->part_id << " of file \'" << el->get_filename() << "\' ";
  oss << "(" << ret/(1UL << 20) << " MB)";
  print_log(oss.str());

  el->is_data = ret > 0;
  
  el->mtx->unlock();
}

} // namespace graphee
