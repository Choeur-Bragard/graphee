#ifndef GRAPHEE_UTILS_H__
#define GRAPHEE_UTILS_H__

#include <iostream>
#include <sstream>
#include <string>
#include <cstring>
#include <ctime>
#include <iomanip>

namespace graphee
{

class utils
{
public:
  enum
  {
    PLAIN,
    GZ,
    BIN,
    SNAPPY
  };
  static const int DIRECT = 0x00000001;
  static const int TRANS = 0x00000010;
  static const int IB = 0x00000100;
  static const int OB = 0x00001000;
};

void print_log(std::string message);
void print_strong_log(std::string message);
void print_warning(std::string message);
void print_error(std::string message);

} // namespace graphee

#endif // GRAPHEE_UTILS_H__
