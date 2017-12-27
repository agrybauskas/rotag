package AtomProperties;

use strict;
use warnings;

use Exporter qw( import );
our @EXPORT_OK = qw( %ATOMS );

# ------------------------------- Atom properties ----------------------------- #

#
# Stores atom properties, parameters, such as, Van der Waals radii, covalent
# radii and bond length.
#

# TODO: should determine covalent radii by investigating amino acid structures.
our %ATOMS = (
    "." => { # Placeholder atom, point in space.
	     "covalent_radius" => {
		                    "length"    => [ 1 ],
				    "error"     => [ 1 ]
				  },
	     "lone_pairs" => 0,
	     "vdw_radius" => 1,
             "valence" => 1,
             "partial_charge" => 0
           },
    "H" => {
	     "covalent_radius" => {
		                    "length"    => [ 0.31 ],
				    "error"     => [ 0.05 ]
				  },
	     "lone_pairs" => 0,
	     "vdw_radius" => 1.2,
             "valence" => 1,
             "partial_charge" => 0.27
           },
    "C" => {
	     "covalent_radius" => {
		                    "length"    => [ 0.76, 0.70, 0.64 ],
				    "error"     => [ 0.03, 0.03, 0.03 ]
				  },
	     "lone_pairs" => 0,
	     "vdw_radius" => 1.77,
             "valence" => 4,
             "partial_charge" => 0.12,
           },
    "N" => {
	     "covalent_radius" => {
		                    "length"    => [ 0.71, 0.60 ],
				    "error"     => [ 0.02, 0.06 ]
				  },
	     "lone_pairs" => 1,
	     "vdw_radius" => 1.66,
             "valence" => 3,
             "partial_charge" => -0.42
           },
    "O" => {
	     "covalent_radius" => {
		                    "length"    => [ 0.66, 0.57 ],
				    "error"     => [ 0.02, 0.07 ]
				  },
	     "lone_pairs" => 2,
	     "vdw_radius" => 1.5,
             "valence" => 2,
             "partial_charge" => -0.57,
           },
    "S" => {
	     "covalent_radius" => {
		                    "length"    => [ 1.05 ],
				    "error"     => [ 0.03 ]
				  },
	     "lone_pairs" => 2,
	     "vdw_radius" => 1.89,
             "valence" => 2,
             "partial_charge" => 0
           },
    "P" => {
	     "covalent_radius" => {
		                    "length"    => [ 1.07 ],
				    "error"     => [ 0.03 ]
				  },
	     "lone_pairs" => 0,
	     "vdw_radius" => 1.9,
             "valence" => 5,
             "partial_charge" => 0
           }
    );

1;
