#include "edgelist.hpp"

namespace graphee
{

bool Edgelist::read_file(std::stringstream &sstream_read)
{
  sstream_read.str("");
  sstream_read.clear();

  // C style code zone, NEEDED FOR ZLIB
  if (file_ptr == nullptr) // not initialized yet
  {
    filelist_it = filelist.begin();
    file_ptr = fopen(filelist_it->c_str(), "rb");
  }
  else if (file_ptr)
  {
    if (feof(file_ptr))
    {
      fclose(file_ptr);
      filelist_it++;
      if (filelist_it != filelist.end())
        file_ptr = fopen(filelist_it->c_str(), "rb");
      else
        return false;
    }
  }
  else
  {
    print_error("Cannot open edgelist file");
    exit(-1);
  }
  // End of C style code zone

  std::thread t (deflate_chunk, this, std::ref(sstream_read));
  t.detach();

  return true;
}


void Edgelist::deflate_chunk(Edgelist* el, std::stringstream &sstream_read)
{
  /** Inspired from : https://zlib.net/zpipe.c
   *
   *  Mixing C and C++, because of the C API of zlib
   */

  int ret;
  size_t have;
  size_t ss_len {0};
  z_stream strm;

  const size_t buf_size {1UL << 18}; // 256kB
  std::array<unsigned char, buf_size> inbuf;
  std::array<unsigned char, buf_size> outbuf;

  /* allocate inflate state */
  strm.zalloc = Z_NULL;
  strm.zfree = Z_NULL;
  strm.opaque = Z_NULL;
  strm.avail_in = 0;
  strm.next_in = Z_NULL;
  ret = inflateInit(&strm);
  if (ret != Z_OK)
    return;

  /** decompress until deflate stream ends or end of file
   * done when inflate() says it's done
   */
  while(ret != Z_STREAM_END && ss_len < el->sstream_limit)
  {
    strm.avail_in = fread(inbuf.data(), 1, buf_size, el->file_ptr);
    if (ferror(el->file_ptr))
    {
      (void)inflateEnd(&strm);
      return;
    }
    if (strm.avail_in == 0)
      break;
    strm.next_in = inbuf.data();

    /* run inflate() on input until output buffer not full */
    while(strm.avail_out == 0)
    {
      strm.avail_out = buf_size;
      strm.next_out = outbuf.data();
      ret = inflate(&strm, Z_NO_FLUSH);
      assert(ret != Z_STREAM_ERROR);  /* state not clobbered */
      switch (ret)
      {
      case Z_NEED_DICT:
        ret = Z_DATA_ERROR;     /* and fall through */
      case Z_DATA_ERROR:
      case Z_MEM_ERROR:
        (void)inflateEnd(&strm);
        return;
      }
      have = buf_size - strm.avail_out;
      sstream_read << ((char*) outbuf.data());
      ss_len += have;
    }
  }

  /* clean up and return */
  (void)inflateEnd(&strm);
  return;
}

} // namespace graphee
