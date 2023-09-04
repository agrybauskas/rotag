#ifndef _PARAMETERS_H_
#define _PARAMETERS_H_

#include <string>

class Parameters {
  public:
    const double EPSILON = epsilon();
    const double PI = pi();

    Parameters(char* parameter_file);
    ~Parameters();

  private:
    double epsilon();
    double pi();
};

#endif
