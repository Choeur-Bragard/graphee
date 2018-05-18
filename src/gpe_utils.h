#ifndef GPE_UTILS_H
#define GPE_UTILS_H

#include <iostream>
#include <string>
#include <cstring>

#include <snappy.h>

bool compress_snappy (const char* in_data, size_t in_bytes, char** out_data, size_t& out_bytes);
bool uncompress_snappy (const char* in_data, size_t in_bytes, char* out_data, size_t& out_bytes);

#endif // GPE_UTILS_H
