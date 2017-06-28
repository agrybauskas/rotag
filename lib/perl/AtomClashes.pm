package AtomClashes;

use strict;
use warnings;

use Exporter qw( import );
our @EXPORT_OK = qw( radius_only );

use List::Util qw( max );

use lib qw( ./ );
use LoadParams qw( covalent_radii
                   vdw_radii );
use Data::Dumper;
my $vdw_file = "../../parameters/vdw_radii.csv";
my $covalent_file = "../../parameters/covalent_radii.csv";

# --------------------------- Detection of atom clashes ----------------------- #

#
# Checks if atoms have clashes with other atoms and removes if they do.
#

#
# Parameters.
#

my %VDW_RADII = %{ vdw_radii( $vdw_file ) };
my $VDW_MAX = max( map { $VDW_RADII{$_} } keys( %VDW_RADII ) ) * 2;

#
# Simplest function for determining atoms clashes. Only radius of atoms are
# considered.
#

sub radius_only
{
    my ( $atom_site, $atom_specifier ) = @_;

    # Clashes of all atoms analyzed, if no specific atoms are selected.
    $atom_specifier = { "group_pdb" => [ "ATOM" ] } unless $atom_specifier;


}

1;
