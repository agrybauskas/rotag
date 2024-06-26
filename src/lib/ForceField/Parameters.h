#ifndef SRC_LIB_FORCEFIELD_PARAMETERS_H_
#define SRC_LIB_FORCEFIELD_PARAMETERS_H_

#include <cmath>
#include <map>
#include <iostream>
#include <string>
#include <vector>

#include <boost/filesystem.hpp>

extern "C" {
    #include "cif.h"
    #include "cif_compiler.h"
    #include "cifvalue.h"
    #include "datablock.h"
}

#include "../PDBx.h"
#include "../Combinatorics.h"

const std::vector<std::string> PARAMETER_TAGS = {
    // "_rotag_parameters" category-related.
    "_rotag_force_field.lj_k",
    "_rotag_force_field.c_k",
    "_rotag_force_field.h_k",
    "_rotag_force_field.t_k",
    "_rotag_force_field.cutoff_atom",
    "_rotag_force_field.cutoff_start",
    "_rotag_force_field.cutoff_end",
    "_rotag_atom_properties.type_symbol",
    "_rotag_atom_properties.hybridization",
    "_rotag_atom_properties.covalent_radius_value",
    "_rotag_atom_properties.covalent_radius_error",
    "_rotag_atom_properties.vdw_radius",
    "_rotag_atom_properties.lone_pair_count",
    "_rotag_atom_properties.valence",
    "_rotag_lennard_jones.type_symbol_1",
    "_rotag_lennard_jones.type_symbol_2",
    "_rotag_lennard_jones.sigma",
    "_rotag_lennard_jones.epsilon",
    "_rotag_partial_charge.label_comp_id",
    "_rotag_partial_charge.label_atom_id",
    "_rotag_partial_charge.value",
    "_rotag_torsional_atom_names.label_comp_id",
    "_rotag_torsional_atom_names.label_atom_id",
    "_rotag_torsional_atom_names.alt_atom_name",
    "_rotag_torsional.label_atom_1_id",
    "_rotag_torsional.label_atom_2_id",
    "_rotag_torsional.label_atom_3_id",
    "_rotag_torsional.label_atom_4_id",
    "_rotag_torsional.epsilon",
    "_rotag_torsional.phase",
    "_rotag_torsional.gamma",
    "_rotag_h_bond.type_symbol",
    "_rotag_h_bond.sigma",
    "_rotag_h_bond.epsilon",
    "_rotag_residue_atom_necessity.label_comp_id",
    "_rotag_residue_atom_necessity.label_atom_id",
    "_rotag_residue_atom_necessity.value",
    "_rotag_clear_hybridization.label_comp_id",
    "_rotag_clear_hybridization.label_atom_id",
    "_rotag_clear_hybridization.type",
    "_rotag_connectivity.label_comp_id",
    "_rotag_connectivity.label_atom_1_id",
    "_rotag_connectivity.label_atom_2_id",
    "_rotag_hydrogen_names.label_comp_id",
    "_rotag_hydrogen_names.label_atom_id",
    "_rotag_hydrogen_names.label_hydrogen_atom_id",
    "_rotag_symmetrical_atom_names.label_comp_id",
    "_rotag_symmetrical_atom_names.label_atom_1_id",
    "_rotag_symmetrical_atom_names.label_atom_2_id",
    "_rotag_dihedral_angle.label_comp_id",
    "_rotag_dihedral_angle.angle",
    "_rotag_dihedral_angle.range_from",
    "_rotag_dihedral_angle.range_to",
    "_rotag_dihedral_angle.step",
    "_rotag_dihedral_angle.type",
    "_rotag_interaction_atom_names.label_atom_id",
    "_rotag_mainchain_atom_names.label_atom_id",
    "_rotag_sidechain_atom_names.label_atom_id",
    "_rotag_rotatable_residue_names.label_comp_id"
};

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
    int64_t lone_pair_count;
    int64_t valence;
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
