// -*- mode: c++; -*-

#ifndef __SIMPLE_EXCEPTION__
#define __SIMPLE_EXCEPTION__

#include <exception>
#include <string>

using std::string;
using std::exception;

// Simple Exception.
class SimpleException : public exception
{
  const string summary;

public:

  SimpleException(const char* desc) : summary(desc) {};
  SimpleException(const string& desc) : summary(desc) {};

  const char* what() const throw()
  {
    return summary.c_str();
  }

};

#endif
