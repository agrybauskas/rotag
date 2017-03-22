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

sub find_flex_bonds # flex - flexible
{
    my %atom_connections = @_;

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
    follow_bond( %atom_connections );

}

1;
