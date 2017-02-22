package SidechainModels;

use strict;
use warnings;

use lib qw( ./ );
use CifParser;
use ConnectAtoms;
use AlterMolecule;

use Data::Dumper;

# ---------------------- Constraint-free Sidechain models --------------------- #

#
# Discribes sidechain models that are not restrained by van der Waals radius.
# These models are idealistic and more developmental than final.
#

#
# Model of serine sidechain.
# Input  (1 arg):
# Output (1 arg):
#

sub cf_serine # cf - abbreviation for constraint-free.
{
    my @amino_acid_data = @_;

    my @atom_connections =
	ConnectAtoms::connect_atoms( 1.51, # HACK: approximate bond-length
				     0.15, # HACK: approximate bond-length-error
				     @amino_acid_data );    
}

1;
