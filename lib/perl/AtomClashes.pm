package AtomClashes;

use strict;
use warnings;

use Exporter qw( import );
our @EXPORT_OK = qw( radius_only );

use lib qw( ./ );

my $parameter_file = "../../parameters/vdw_radii.csv";

# --------------------------- Detection of atom clashes ----------------------- #

#
# Checks if atoms have clashes with other atoms and removes if they do.
#

#
# Parameters.
#

#
# Converts parameter file that identifies rotatable side-chain bonds to hash of
# hashes. Ex.:
# {
#   "N" => 1.66,
#   "C" => 1.77
# }
#

my %VDW_RADII;
my $atom_name;
my $atom_radius;

{
    open( my $fh, "<", $parameter_file )
    	or die "Can't open < van_der_waals_radii.csv: $!";
    foreach( <$fh> ) {
	chomp( $_ );
	( $atom_name, $atom_radius ) = split( ",", );
	$VDW_RADII{"$atom_name"} = $atom_radius;
    }
    close($fh);
}

#
# Simplest function for determining atoms clashes. Only radius of atoms are
# considered.
#

sub radius_only
{

}
