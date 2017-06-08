package LoadParams;

use strict;
use warnings;

use Exporter qw( import );
our @EXPORT_OK = qw( rotatable_bonds
                     vdw_radius );

#
# Converts parameter file that identifies rotatable side-chain bonds to hash of
# hashes.
# Input  (1 arg): destination of csv parameter file.
# Output (1 arg): data structure: hashes of hashes of multidimentional array.
# Ex.:
# {
#   "SER" => { OG => [ [ "CA", "CB" ] ] }
# }
#

sub rotatable_bonds
{
    my ( $parameter_file ) = @_;

    my %ROTATABLE_BONDS;

    open( my $fh, "<", $parameter_file )
    	or die "Can't open < rotatable_bonds.csv: $!";

    for my $data_row ( map { [ split( ",", $_ ) ] } <$fh> ) {
	$ROTATABLE_BONDS{$data_row->[0]}{$data_row->[1]} = [];
	for my $bond ( @{ $data_row }[2..$#{ $data_row }] ) {
	    push( @{ $ROTATABLE_BONDS{$data_row->[0]}{$data_row->[1]} },
		  [ split( ":", $bond ) ] );
	}
    }

    return \%ROTATABLE_BONDS;
}

#
# Converts parameter file that identifies atom Van der Waals radius to simple
# hash.
# Input  (1 arg): destination of csv parameter file.
# Output (1 arg): hash of values of atom vdw radii.
# Ex.:
# {
#   "N" => 1.66,
#   "C" => 1.77
# }
#

sub vdw_radius
{
    my $parameter_file;

    my %VDW_RADII;
    my $atom_name;
    my $atom_radius;

    open( my $fh, "<", $parameter_file )
    	or die "Can't open < van_der_waals_radii.csv: $!";
    foreach( <$fh> ) {
	chomp( $_ );
	( $atom_name, $atom_radius ) = split( ",", );
	$VDW_RADII{"$atom_name"} = $atom_radius;
    }
    close($fh);

    return \%VDW_RADII;
}

1;
