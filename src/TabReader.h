/*******************************************************************************
 Copyright 2006-2012 Lukas Käll <lukas.kall@scilifelab.se>

 Licensed under the Apache License, Version 2.0 (the "License");
 you may not use this file except in compliance with the License.
 You may obtain a copy of the License at

 http://www.apache.org/licenses/LICENSE-2.0

 Unless required by applicable law or agreed to in writing, software
 distributed under the License is distributed on an "AS IS" BASIS,
 WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 See the License for the specific language governing permissions and
 limitations under the License.

 *******************************************************************************/
#ifndef TABREADER_H_
#define TABREADER_H_

#include <string>
#include <string_view>
#include <charconv>
#include <climits>
#include <cstring>
#include <cmath>

// using char pointers is much faster than istringstream
class TabReader {
 public:
  TabReader(const std::string& line)
      : f_(line.c_str()), end_(line.c_str() + line.size()), err(0) {
    errno = 0;
  }
  
  void advance(const char* next) {
    if (*next != '\0') {
      f_ = next + 1; // eats up the tab
    } else {
      f_ = next; // prevents pointing over the null byte
    }
  }
  
  void skip() {
    const char* pch = strchr(f_, '\t');
    if (pch == NULL) {
      err = 1;
    } else {
      advance(pch);
    }
  }
  
  void skip(size_t numSkip) {
    for (size_t i = 0; i < numSkip; ++i) skip();
  }
  
  double readDouble() {
#if defined(__cpp_lib_to_chars)
    // std::from_chars rounds correctly like strtod but is much faster. Any
    // input it does not fully accept (leading '+' or whitespace, hex,
    // out-of-range values, trailing garbage) takes the strtod path below, so
    // results and error handling are unchanged.
    {
      double d;
      std::from_chars_result res = std::from_chars(f_, end_, d);
      if (res.ec == std::errc() && res.ptr != f_ &&
          (*res.ptr == '\0' || isspace(static_cast<unsigned char>(*res.ptr)))) {
        advance(res.ptr);
        return d;
      }
    }
#endif
    char* next = NULL;
    errno = 0;
    double d = strtod(f_, &next);
    if (next == f_ || (*next != '\0' && !isspace(static_cast<unsigned char>(*next)))
                   || ((d == HUGE_VAL || d == -HUGE_VAL) && errno == ERANGE)) {
      err = errno ? errno : 1;
    }
    advance(next);
    return d;
  }
  
  int readInt() {
    char* next = NULL;
    errno=0;
    long val = strtol(f_, &next, 10);
    if (next == f_ || (*next != '\0' && !isspace(static_cast<unsigned char>(*next)))
                   || val < INT_MIN || val > INT_MAX) {
      err = errno ? errno : 1;
    }
    advance(next);
    return static_cast<int>(val);
  }
  
  std::string readString() {
    return std::string(readView());
  }

  // the next field, as a view into the line
  std::string_view readView() {
    const char* pch = strchr(f_, '\t');
    if (pch == NULL) {
      err = 1;
      return std::string_view(f_);
    } else {
      std::string_view s(f_, static_cast<size_t>(pch - f_));
      advance(pch);
      return s;
    }
  }

  bool error() { return err != 0; }
 private:
  const char* f_;
  const char* end_;
  int err;
};

#endif /*TABREADER_H_*/