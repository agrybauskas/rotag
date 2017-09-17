package AtomProperties;

use strict;
use warnings;

use Exporter qw( import );
our @EXPORT_OK = qw( ATOMS );

# ------------------------------- Atom properties ----------------------------- #

#
# Stores atom properties, parameters, such as, Van der Waals radii, covalent
# radii and bond length.
#

our %ATOMS = (
    "H" => {
	     "covalent_radius" => {
		                    "length" => [ 0.31 ],
				    "error"  => [ 0.15 ]
				  },
	     "vdw_radius" => 1.2
           },
    "C" => {
	     "covalent_radius" => {
		                    "length" => [ 0.76, 0.73, 0.69 ],
				    "error"  => [ 0.03, 0.06, 0.03 ]
				  },
	     "vdw_radius" => 1.77
           },
    "N" => {
	     "covalent_radius" => {
		                    "length" => [ 0.71, 0.6 ],
				    "error"  => [ 0.03, 0.03 ]
				  },
	     "vdw_radius" => 1.66
           },
    "O" => {
	     "covalent_radius" => {
		                    "length" => [ 0.66, 0.57 ],
				    "error"  => [ 0.06, 0.06 ]
				  },
	     "vdw_radius" => 1.5
           },
    "S" => {
	     "covalent_radius" => {
		                    "length" => [ 1.05 ],
				    "error"  => [ 0.09 ]
				  },
	     "vdw_radius" => 1.89
           },
    "P" => {
	     "covalent_radius" => {
		                    "length" => [ 1.07 ],
				    "error"  => [ 0.09 ]
				  },
	     "vdw_radius" => 1.9
           }
    )

1;
