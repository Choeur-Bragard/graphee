#include "gpe_utils.h"

void compress_snappy (char* in_data, size_t in_bytes, char* out_data, size_t& out_bytes) {
  char* uncomp_data = in_data;
  size_t max_bytes_per_block {UINT32_MAX - 1};

  int nblocks = in_bytes/max_bytes_per_block + (in_bytes%max_bytes_per_block == 0 ? 0 : 1);

  size_t max_size = nblocks*sizeof(size_t) + nblocks*snappy::MaxCompressedLength(max_bytes_per_block);

  char* comp_data = out_data; 

  std::memcpy (comp_data, &nblocks, sizeof(nblocks));
  comp_data += nblocks;

  bool comp_succeed {false};
  size_t read_bytes {0};
  size_t comp_size;
  size_t uncomp_size;

  out_bytes = 0;

  for (int blkID = 0; blkID < nblocks; blkID++) {
    uncomp_size = std::min(max_bytes_per_block, in_bytes-read_bytes);

    snappy::RawCompress (uncomp_data, uncomp_size, comp_data, &comp_size);

    read_bytes += uncomp_size;
    uncomp_data += uncomp_size;

    comp_data += comp_size;
    out_bytes += comp_size;

    if (uncomp_data > in_data + in_bytes || comp_data > out_data + max_size) {
      gpe_error ("Snappy compression failed");
      exit (-1);
    }
  }
}

bool uncompress_snappy (char* in_data, size_t in_bytes, char* out_data, size_t out_bytes) {
  char* comp_data = in_data;
  char* uncomp_data = out_data;

  int nblocks;
  std::memcpy (&nblocks, comp_data, sizeof(nblocks));

  comp_data += sizeof(nblocks);

  size_t comp_size;
  size_t uncomp_size;

  bool uncomp_succeed {false};

  for (int blkID = 0; blkID < nblocks; blkID++) {
    std::memcpy (&comp_size, comp_data, sizeof(comp_size));
    comp_data += sizeof(comp_data);

    snappy::GetUncompressedLength (comp_data, comp_size, &uncomp_size);
    uncomp_succeed = snappy::RawUncompress (comp_data, comp_size, uncomp_data);
    if (!uncomp_succeed) {
      return false;
    }

    comp_data += comp_size;
    uncomp_data += uncomp_size;

    if (comp_data > in_data+in_bytes || uncomp_data > out_data+out_bytes) {
      return false;
    }
  }

  return true;
}

void gpe_log (std::string message) {
  std::time_t t = std::time(nullptr);
  std::tm tm = *std::localtime(&t);

  std::cout << "\033[1m [" << std::put_time(&tm, "%T") << "] [GRAPHEE] [LOG]  \033[0m" << message << std::endl;
}

void gpe_warning (std::string message) {
  std::time_t t = std::time(nullptr);
  std::tm tm = *std::localtime(&t);

  std::cout << "\033[1;36m [" << std::put_time(&tm, "%T") << "] [GRAPHEE] [WARN]\033[0m \033[36m" << message << "\033[0m" << std::endl;
}

void gpe_error (std::string message) {
  std::time_t t = std::time(nullptr);
  std::tm tm = *std::localtime(&t);

  std::cout << "\033[1;31m [" << std::put_time(&tm, "%T") << "] [GRAPHEE] [ERR] \033[0m \033[31m" << message << "\033[0m" << std::endl;
}
