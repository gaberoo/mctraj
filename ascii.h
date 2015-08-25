#ifndef __ASCII_H__
#define __ASCII_H__

#include <string>

namespace ascii {
  const string red     = "\033[00;31m";
  const string green   = "\033[00;32m";
  const string yellow  = "\033[00;33m";
  const string blue    = "\033[00;34m";
  const string magenta = "\033[00;35m";
  const string cyan    = "\033[00;36m";
  const string white   = "\033[00;37m";

  const string b_red     = "\033[01;31m";
  const string b_green   = "\033[01;32m";
  const string b_yellow  = "\033[01;33m";
  const string b_blue    = "\033[01;34m";
  const string b_magenta = "\033[01;35m";
  const string b_cyan    = "\033[01;36m";
  const string b_white   = "\033[01;37m";

  const string end     = "\033[00m";
}

#endif
