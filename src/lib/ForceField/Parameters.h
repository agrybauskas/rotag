#ifndef _PARAMETERS_H_
#define _PARAMETERS_H_

#include <map>
#include <string>

struct CovalentRadius {
  double value;
  double error;
};

struct LennardJones {
  double sigma;
  double epsilon;
};

struct PartialCharge {
  double value;
};

struct AtomProperties {
  std::map<std::string, CovalentRadius> covalent_radius;
  double vdw_radius;
  double lone_pair_count;
  double valence;
};

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

    double lj_k;
    double c_k;
    double h_k;
    double t_k;
    double cutoff_atom;
    double cutoff_start;
    double cutoff_end;

    std::map<std::string, AtomProperties> ATOM_PROPERTIES;
    std::map<std::string, std::map<std::string, LennardJones>> LENNARD_JONES;
    std::map<std::string, std::map<std::string, PartialCharge>> PARTIAL_CHARGE;

    double max_connection_length;
    double max_interaction_length;

    Parameters(char* program_file_path);
    ~Parameters();

  private:
    double epsilon();
    double pi();
};

#endif
