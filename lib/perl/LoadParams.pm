package LoadParams;

use strict;
use warnings;

use Exporter qw( import );
our @EXPORT_OK = qw( covalent_radii rotatable_bonds vdw_radii );

#
# Converts parameter file that identifies rotatable side-chain bonds to special
# data structure.
# Input:
#     $parameter_file - destination of csv parameter file describing rotatable
#     bonds in specific side-chains.
# Output:
#     %ROTATABLE_BONDS - special data structure describing rotatable bonds of
#     each side-chain.
#     Ex.: { "SER" => { OG => [ [ "CA", "CB" ] ] } }
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
# Input:
#     $parameter_file - destination of csv parameter file describing vdw radius
#     of specific atom.
# Output:
#     %VDW_radii - hash describing atom vdw radii.
#     Ex.: {
#            "N" => 1.66,
#            "C" => 1.77
#          }
#

sub vdw_radii
{
    my ( $parameter_file ) = @_;

    my %VDW_RADII;
    my $atom_name;
    my $atom_radius;

    open( my $fh, "<", $parameter_file )
    	or die "Can't open < vdw_radii.csv: $!";
    for my $row ( <$fh> ) {
    	chomp( $row );
    	( $atom_name, $atom_radius ) = split( ",", $row );
    	$VDW_RADII{"$atom_name"} = $atom_radius;
    }
    close($fh);

    return \%VDW_RADII;
}

#
# Converts parameter file that identifies bond length and converts to a specific
# data structure.
# Input:
#     $parameter_file - destination of csv parameter file describing covalent
#     radii and its error.
# Output:
#     %COVALENT_RADII - special data structure describing atom bond lengths and
#     bond length errors.
#     Ex.: {
#            "C" => "length" => [ 1.513, 1.455 ],
#                   "error"  => [ 0.014, 0.011 ]
#          }
#

sub covalent_radii
{
    my ( $parameter_file ) = @_;

    my %COVALENT_RADII;
    my $atom_name;
    my $bond_length;
    my $length_error;

    open( my $fh, "<", $parameter_file )
    	or die "Can't open < vdw_radii.csv: $!";
    foreach( <$fh> ) {
    	chomp( $_ );
	( $atom_name, $bond_length, $length_error ) =
	    split( ",", $_ );
	if( exists $COVALENT_RADII{$atom_name} ) {
	    push( @{ $COVALENT_RADII{$atom_name}{"bond_length"} },
		  $bond_length );
	    push( @{ $COVALENT_RADII{$atom_name}{"length_error"} },
		  $length_error );
	} else {
	    $COVALENT_RADII{$atom_name}{"bond_length"} = [ $bond_length ];
	    $COVALENT_RADII{$atom_name}{"length_error"} = [ $length_error ];
	}
    }
    close($fh);

    return \%COVALENT_RADII;
}

1;
