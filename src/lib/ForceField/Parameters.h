#ifndef SRC_LIB_FORCEFIELD_PARAMETERS_H_
#define SRC_LIB_FORCEFIELD_PARAMETERS_H_

#include <map>
#include <string>
#include <vector>

struct AltName {
  std::string alt_name;
};

struct BondType {
  double min_length;
  double max_length;
};

struct ClearHybridization {
  std::string type;
};

struct CovalentBondCombinations {
  std::vector<std::vector<double>> values;
  std::vector<std::vector<double>> errors;
};

struct CovalentRadius {
  double value;
  double error;
};

struct DihedralAngleRestraint {
  double range_from;
  double range_to;
  double step;
  std::string type;
};

struct HBond {
  double epsilon;
  double phase;
  double gamma;
};

struct LennardJones {
  double sigma;
  double epsilon;
};

struct PartialCharge {
  double value;
};

struct Torsional {
  double epsilon;
  double phase;
  double gamma;
};

struct AtomProperties {
  std::map<std::string, CovalentRadius> covalent_radius;
  double vdw_radius;
  int lone_pair_count;
  int valence;
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
    std::map<std::string, std::map<std::string, AltName>> TORSIONAL_ATOM_NAMES;
    std::map<std::string, Torsional> TORSIONAL;
    std::map<std::string, HBond> H_BOND;
    std::map<std::string, std::map<std::string, bool>> RESIDUE_ATOM_NECESSITY;
    std::map<std::string, std::map<std::string,
                                   ClearHybridization>> CLEAR_HYBRIDIZATION;
    std::map<std::string, std::map<std::string,
                                   std::vector<std::string>>> CONNECTIVITY;
    std::map<std::string, std::map<std::string,
                                   std::vector<std::string>>> HYDROGEN_NAMES;
    std::map<std::string,
             std::map<std::string,
                      std::vector<std::string>>> SYMMETRICAL_ATOM_NAMES;
    std::map<std::string, std::map<std::string,
                                   DihedralAngleRestraint>> DIHEDRAL_ANGLE;
    std::vector<std::string> INTERACTION_ATOM_NAMES;
    std::vector<std::string> MAINCHAIN_ATOM_NAMES;
    std::vector<std::string> SIDECHAIN_ATOM_NAMES;
    std::vector<std::string> ROTATABLE_RESIDUE_NAMES;
    std::map<std::string, std::map<std::string, std::map<std::string,
                                                         BondType>>> BOND_TYPE;
    std::map<std::string, std::vector<double>> COVALENT_RADII_VALUES;
    std::map<std::string, std::vector<double>> COVALENT_RADII_ERRORS;
    std::map<std::string,
             std::map<std::string,
                      CovalentBondCombinations>> COVALENT_BOND_COMBINATIONS;

    double max_connection_length = 0.0;
    double max_interaction_length = 0.0;

    explicit Parameters(char* program_file_path);
    ~Parameters();

 private:
    double epsilon();
    double pi();
};

#endif  // SRC_LIB_FORCEFIELD_PARAMETERS_H_
