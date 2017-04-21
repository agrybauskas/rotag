package SidechainModels;

use strict;
use warnings;

use lib qw( ./ );
use CifParser;
use ConnectAtoms;
use AlterMolecule;
use LinearAlgebra;

use feature qw( current_sub );
use Data::Dumper;

my $parameter_file = "../../parameters/rotatable_bonds.csv";

# ------------------------ Idealistic sidechain models ------------------------ #

#
# Parameters
#

#
# Converts parameter file that identifies rotatable side-chain bonds to hash of
# hashes. Ex.:
# {
#   "SER" => { OG => [ [ "CA", "CB" ] ] }
# }
#

my %ROTATABLE_BONDS;

{
    open( my $fh, "<", $parameter_file )
    	or die "Can't open < rotatable_bonds.csv: $!";

    for my $data_row ( map { [ split( ",", $_ ) ] } <$fh> ) {
	$ROTATABLE_BONDS{$data_row->[0]}{$data_row->[1]} = [];
	for my $bond ( @{ $data_row }[2..$#{ $data_row }] ) {
	    push( @{ $ROTATABLE_BONDS{$data_row->[0]}{$data_row->[1]} },
		  [ split( ":", $bond ) ] );
	}
    }
}

#
# Discribes sidechain models that are not restrained by van der Waals radius.
# These models are idealistic and more developmental than final. One model per
# residue.
#

#
# Model that uses only rotation around single bonds.
# Input  (1 arg):
# Output (1 arg):
#

sub rotation_only
{
    my ( $atom_site, $atom_specifier ) = @_;

    # Selects specified atom(s) id(s).
    my $target_atom_id =
	&CifParser::select_atom_data( [ "id" ],
				      &CifParser::filter_atoms( $atom_specifier,
								$atom_site ) );

    # Iterates through target atom(s) and assigns conformational equations which
    # can produce pseudo-atoms later.
    my @rotatable_bonds;

    for my $id ( @$target_atom_id ) {
	@rotatable_bonds = @{ $ROTATABLE_BONDS{"SER"}{"OG"} };
	# Creates matrices for atom alterations.
	for( my $i = 0; $i < scalar( @rotatable_bonds ); $i++ ) {
	}
    }
}

1;
