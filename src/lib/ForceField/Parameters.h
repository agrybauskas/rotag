#ifndef _PARAMETERS_H_
#define _PARAMETERS_H_

class Parameters {
  public:
    const double EPSILON = epsilon();
    const double PI = pi();
    const char* SIG_FIG_MIN = "%.3f";
    const char* SIG_FIG_MAX = "%.6f";

    Parameters(char* parameter_file);
    ~Parameters();

  private:
    double epsilon();
    double pi();
};

#endif
