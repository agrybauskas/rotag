package SidechainModels;

use strict;
use warnings;

use lib qw( ./ );
use CifParser;
use ConnectAtoms;
use AlterMolecule;

use Data::Dumper;

my $parameter_file = "../../parameters/rotatable_bonds.tsv";

# ------------------------ Idealistic sidechain models ------------------------ #

#
# Parameters
#

#
# Converts parameter file that identifies rotatable side-chain bonds to array of
# arrays. Ex. ( [ "SER", "N",  "CA", "no" ],
#               [ "SER", "CA", "C",  "no" ],
#               [ "SER", "C",  "O",  "no" ],
#               [ "SER", "CA", "CB", "yes"],
#               [ "SER", "CB", "OG", "no" ] )
#

my @ROTATABLE_BONDS;

{
    open( my $fh, "<", $parameter_file )
	or die "Can't open < rotatable_bonds.tsv: $!";
    # TODO: instead of list of lists, make hash of lists. Might be more effective
    #       method to search.
    foreach( map { [ split(" ", $_) ] }
	     grep { /^\w+\s+\w+\s+\w+\s+(yes|no)$/ } <$fh> ) {
	push( @ROTATABLE_BONDS, $_ )
    }
    close( $fh );
}

#
# Discribes sidechain models that are not restrained by van der Waals radius.
# These models are idealistic and more developmental than final.
#

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
    my $atom_site = shift;

    # Connects atoms.
    my $connected_atoms =
	ConnectAtoms::connect_atoms( 1.592, # HACK: empirical bond length.
				     0.404, # HACK: empirical bond length error.
				     $atom_site );
}

1;
