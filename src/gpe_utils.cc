#include "gpe_utils.h"

size_t max_compress_size (size_t in_bytes) {
  const size_t max_bytes_per_block {UINT32_MAX - 1};
  int nblocks = in_bytes/max_bytes_per_block + (in_bytes%max_bytes_per_block == 0 ? 0 : 1);
  return snappy::MaxCompressedLength(in_bytes) + nblocks*sizeof(size_t);
}

void compress_snappy (char* in_data, size_t in_bytes, char* out_data, size_t& out_bytes) {
  const size_t max_bytes_per_block {UINT32_MAX - 1};

  size_t in_pos {0};
  size_t out_pos {0};

  int nblocks = in_bytes/max_bytes_per_block + (in_bytes%max_bytes_per_block == 0 ? 0 : 1);

  std::memcpy (out_data + out_pos, &nblocks, sizeof(nblocks));
  out_pos += sizeof(nblocks);

  size_t in_block_size;
  size_t out_block_size;

  for (int blockID = 0; blockID < nblocks; blockID++) {
    in_block_size = std::min(max_bytes_per_block, in_bytes - in_pos);

    snappy::RawCompress (in_data + in_pos, in_block_size, out_data + out_pos + sizeof(size_t), &out_block_size);
    std::memcpy (out_data + out_pos, &out_block_size, sizeof(size_t));

    in_pos += in_block_size;
    out_pos += out_block_size + sizeof(size_t);

    if (in_pos > in_bytes || out_pos > out_bytes) {
      gpe_error ("Snappy compression failed");
      exit (-1);
    }
  }

  out_bytes = out_pos;

  std::ostringstream oss;
  oss << "Snappy compressed " << in_bytes << " into " << out_bytes;
  gpe_log (oss.str());
}

bool uncompress_snappy (char* in_data, size_t in_bytes, char* out_data, size_t out_bytes) {
  size_t in_pos {0};
  size_t out_pos {0};

  int nblocks;

  std::memcpy (&nblocks, in_data, sizeof(nblocks));
  in_pos += sizeof(nblocks);

  std::cout << nblocks << std::endl;

  size_t in_block_size;
  size_t out_block_size;

  bool uncomp_succeed {false};

  for (int blockID = 0; blockID < nblocks; blockID++) {
    std::memcpy (&in_block_size, in_data + in_pos, sizeof(size_t));
    in_pos += sizeof(size_t);
    std::cout << in_block_size << std::endl;

    snappy::GetUncompressedLength (in_data + in_pos, in_block_size, &out_block_size);
    uncomp_succeed = snappy::RawUncompress (in_data + in_pos, in_block_size, out_data + out_pos);
    if (!uncomp_succeed) {
      return false;
    }

    in_pos += in_block_size;
    out_pos += out_block_size;

    if (in_pos > in_bytes || out_pos > out_bytes) {
      std::cout << in_pos << " " << in_bytes << " " << out_block_size << " " << out_bytes << std::endl;
      return false;
    }
  }

  return true;
}

void gpe_log (std::string message) {
  std::time_t t = std::time(nullptr);
  std::tm tm = *std::localtime(&t);

  std::cout << "\033[1m [" << std::put_time(&tm, "%T") << "] [GRAPHEE] [LOG] -\033[0m " << message << std::endl;
}

void gpe_warning (std::string message) {
  std::time_t t = std::time(nullptr);
  std::tm tm = *std::localtime(&t);

  std::cout << "\033[1;33m [" << std::put_time(&tm, "%T") << "] [GRAPHEE] [WRN] -\033[0m \033[93m" << message << "\033[0m" << std::endl;
}

void gpe_error (std::string message) {
  std::time_t t = std::time(nullptr);
  std::tm tm = *std::localtime(&t);

  std::cout << "\033[1;31m [" << std::put_time(&tm, "%T") << "] [GRAPHEE] [ERR] -\033[0m \033[91m" << message << "\033[0m" << std::endl;
}
