#include "gpe_utils.h"

bool compress_snappy (const char* in_data, size_t in_bytes, char** out_data, size_t& out_bytes) {
  char* uncomp_data = in_data;
  size_t max_bytes_per_block {UINT32_MAX - 1};

  int nblocks = in_bytes/max_bytes_per_block + (in_bytes%max_bytes_per_block == 0 ? 0 : 1);

  size_t max_size = nblocks*sizeof(size_t) + nblocks*snappy::MaxCompressedLength(max_bytes_per_block);

  (*out_data) = new char [max_size];
  char* comp_data = (*out_data); 

  std::memcpy (comp_data, nblocks, sizeof(nblocks));
  comp_data += nblocks;

  bool comp_succeed {false};
  size_t read_bytes {0};
  size_t comp_size;
  size_t uncomp_size;

  out_bytes = 0;

  for (int blkID = 0; blkID < nblocks; blkID++) {
    uncomp_size = std::min(max_bytes_per_block, in_bytes-read_bytes);

    comp_succeed = snappy::RawCompress (uncomp_data, uncomp_size, comp_data, &comp_size);
    if (!comp_succeed) {
      return false;
    }

    read_bytes += uncomp_size;
    uncomp_data += uncomp_size;

    comp_data += comp_size;
    out_bytes += comp_size;

    if (uncomp_data > in_data + in_bytes || comp_data > out_data + max_size) {
      return false;
    }
  }
  return true;
}

bool uncompress_snappy (const char* in_data, size_t in_bytes, char* out_data, size_t& out_bytes) {
  char* comp_data = in_data;
  char* uncomp_data = out_data;

  int nblocks;
  std::memcpy (&nblocks, data, sizeof(nblocks));

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

    if (comp_data > in_data+in_bytes || uncomp_data > out_data + out_data) {
      return false;
    }
  }

  return true;
}
