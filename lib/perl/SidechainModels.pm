package SidechainModels;

use strict;
use warnings;

use lib qw( ./ );
use CifParser;
use ConnectAtoms;
use AlterMolecule;

use Data::Dumper;

# ------------------------ Idealistic sidechain models ------------------------ #

#
# Discribes sidechain models that are not restrained by van der Waals radius.
# These models are idealistic and more developmental than final.
#

#
# Recursive function that follows the connections of bonds and returns pairs of
# atom pairs.
# Input   (1 arg):
# Output  (1 arg):
#

sub find_flex_bonds # flex - flexible.
{
    my %atom_connections = @_;

    my @visited_atoms;
    my @flexible_bonds;

    # First, deletes atom connections where atom is connected to only one
    # counterpart. There cannot be meaningful rotation around this type of bond.
    # Also, algorithm always starts from CA. There must be direction of bonds,
    # because matrix product is not commutative.
    for my $atom ( keys %atom_connections ) {
	if( scalar( @{ $atom_connections{$atom} } ) == 1 ) {
	}
    }
}

#
# Model that uses only rotation around single bonds.
# Input  (1 arg): array of arrays containing amino acid data, such as id,
#                 type_symbol, label_alt_id, label_comp_id, Cartn_x, Cartn_y,
#                 Cartn_z:
#                 Ex.: ( [ 151, "CB", "SER", -54.506, -111.705, -5.149 ],
#                        [ 152, "OG", "SER", -54.962, -111.220, -3.896 ] )
# Output (1 arg):
#

sub rotation_only
{
    my @amino_acid_data = @_;

    my %atom_connections =
	%{ ConnectAtoms::connect_atoms( 1.51, # HACK: approximate bond-length
					0.15, # HACK: bond-length-error
					@amino_acid_data ) };

    # Picks any atom in atom_connections hash and start to go through
    # connections.
    find_flex_bonds( %atom_connections );
}

1;
