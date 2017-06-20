package AtomClashes;

use strict;
use warnings;

use Exporter qw( import );
our @EXPORT_OK = qw( radius_only );

use lib qw( ./ );
use LoadParams qw( vdw_radius );

my $parameter_file = "../../parameters/vdw_radii.csv";

# --------------------------- Detection of atom clashes ----------------------- #

#
# Checks if atoms have clashes with other atoms and removes if they do.
#

#
# Parameters.
#

my %VDW_RADII = %{ vdw_radius( $parameter_file ) };

#
# Simplest function for determining atoms clashes. Only radius of atoms are
# considered.
#

# sub radius_only
# {
#     my $atom_site
# }
