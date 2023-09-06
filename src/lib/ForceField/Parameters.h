#ifndef _PARAMETERS_H_
#define _PARAMETERS_H_

class Parameters {
  public:
    const double EPSILON = epsilon();
    const double PI = pi();
    const char* SIG_FIG_MIN = "%.3f";
    const char* SIG_FIG_MAX = "%.6f";
    const double SP3_ANGLE = 109.5 * this->PI / 180.0;
    const double SP2_ANGLE = 120.0 * this->PI / 180.0;
    const double SP_ANGLE = this->PI / 180.0;
    const double SP3_ANGLE_ERR = 5.0 * this->PI / 180.0;
    const double SP2_ANGLE_ERR = 5.0 * this->PI / 180.0;
    const double SP_ANGLE_ERR = 5.0 * this->PI / 180.0;

    Parameters(char* program_file_path);
    ~Parameters();

  private:
    double epsilon();
    double pi();
};

#endif
