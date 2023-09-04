#ifndef _PARAMETERS_H_
#define _PARAMETERS_H_

#include <string>

class Parameters {
  public:
    double PI = pi();

    Parameters(char* parameter_file);
    ~Parameters();

  private:
    double pi();
};

#endif
