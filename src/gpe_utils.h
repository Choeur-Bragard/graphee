#ifndef GPE_UTILS_H
#define GPE_UTILS_H

#include <iostream>
#include <string>
#include <cstring>

#include <snappy.h>

namespace graphee {

class gpe_utils {
public:
  enum {PLAIN, GZ};
  static const int DIRECT = 0x00000001; 
  static const int TRANS  = 0x00000010; 
  static const int UO     = 0x00000100; 
};
}

bool compress_snappy (char* in_data, size_t in_bytes, char** out_data, size_t& out_bytes);
bool uncompress_snappy (char* in_data, size_t in_bytes, char* out_data, size_t& out_bytes);


#endif // GPE_UTILS_H
