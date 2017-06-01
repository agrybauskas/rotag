package AtomClashes;

use strict;
use warnings;

use Exporter qw( import );
our @EXPORT_OK = qw( radius_only );

use lib qw( ./ );

my $parameter_file = "../../parameters/van_der_waals_radii.csv";

# --------------------------- Detection of atom clashes ----------------------- #

#
# Parameters.
#

#
# Converts parameter file that identifies rotatable side-chain bonds to hash of
# hashes. Ex.:
# {
#   "SER" => { OG => [ [ "CA", "CB" ] ] }
# }
#

my %WAALS_RADII;
my $atom_name;
my $atom_radius;

{
    open( my $fh, "<", $parameter_file )
    	or die "Can't open < van_der_waals_radii.csv: $!";
    foreach( <$fh> ) {
	chomp( $_ );
	( $atom_name, $atom_radius ) = split( ",", );
	$WAALS_RADII{"$atom_name"} = $atom_radius;
    }
    close($fh);
}

#
# Checks if atoms have clashes with other atoms and removes if they do.
#

#
# Simplest function for determining atoms clashes. Only radius of atoms are
# considered.
#

sub radius_only
{

}
