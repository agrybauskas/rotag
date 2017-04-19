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
# Converts parameter file that identifies rotatable side-chain bonds to hash of
# hashes.
# Ex. residue  rotatable? 1st atom  2nd atom
# {
#     "SER" => {
# 	         "yes" => {
# 	                    "CA" => [ "CB" ],
# 	                    "CB" => [ "CA" ]
# 	       },
# 	         "no"  => {
# 		            "N"  => [ "CA" ],
# 		            "CA" => [ "N", "C" ],
# 		            "C"  => [ "CA", "O" ],
# 		            "O"  => [ "C" ],
# 		            "CB" => [ "OG" ],
# 		            "OG" => [ "CB" ]
# 	                  }
#              }
# }
#

my %ROTATABLE_BONDS;

{
    open( my $fh, "<", $parameter_file )
    	or die "Can't open < rotatable_bonds.tsv: $!";

    # Checks, if first or second atom is as hash key. If not, adds key and list
    # of connected atoms. If yes, pushes to list of connected atoms.
    foreach( map { [ split(" ", $_) ] }
    	     grep { /^\w+\s+\w+\s+\w+\s+(yes|no)$/ } <$fh> ) {
	# Checks, if first atom is in the hash as key.
	if( ! exists $ROTATABLE_BONDS{$_->[0]}{$_->[3]}{$_->[1]} ) {
	    $ROTATABLE_BONDS{$_->[0]}{$_->[3]}{$_->[1]} = [ $_->[2] ];
	} else {
	    if( exists $ROTATABLE_BONDS{$_->[0]}{$_->[3]}{$_->[1]} ) {
	        push( @{ $ROTATABLE_BONDS{$_->[0]}{$_->[3]}{$_->[1]} },
		      $_->[2] );
	    }
	    if( exists $ROTATABLE_BONDS{$_->[0]}{$_->[3]}{$_->[2]} ) {
		push( @{ $ROTATABLE_BONDS{$_->[0]}{$_->[3]}{$_->[2]} },
		      $_->[1] );
	    }
	}

	# Checks, if second atom is in the hash as key.
	if ( ! exists $ROTATABLE_BONDS{$_->[0]}{$_->[3]}{$_->[2]} ) {
	    $ROTATABLE_BONDS{$_->[0]}{$_->[3]}{$_->[2]} = [ $_->[1] ];
	} else {
	    if( exists $ROTATABLE_BONDS{$_->[0]}{$_->[3]}{$_->[1]} ) {
	        push( @{ $ROTATABLE_BONDS{$_->[0]}{$_->[3]}{$_->[1]} },
		      $_->[2] );
	    }
	    if( exists $ROTATABLE_BONDS{$_->[0]}{$_->[3]}{$_->[2]} ) {
		push( @{ $ROTATABLE_BONDS{$_->[0]}{$_->[3]}{$_->[2]} },
		      $_->[1] );
	    }
	}
    }
}

#
# Discribes sidechain models that are not restrained by van der Waals radius.
# These models are idealistic and more developmental than final.
#

#
# Model that uses only rotation around single bonds.
# Input  (1 arg):
# Output (1 arg):
#

sub rotation_only
{
    my ( $atom_site, $atom_specifier ) = @_;

    # Connects atoms.
    my $connected_atoms =
    	&ConnectAtoms::connect_atoms( 1.592, # HACK: empirical bond length.
				      0.404, # HACK: empirical bond length error.
				      $atom_site );

    # Selects specified atom(s) id(s).
    my $target_atom_id =
	&CifParser::select_atom_data( [ "id" ],
				      &CifParser::filter_atoms(
					  $atom_specifier,
					  $connected_atoms ) );

    # Iterates through target atom(s) and assigns conformational equations which
    # can produce pseudo-atoms later.
    my @rotatable_path; # All rotatable bonds from target atom to CA.

    for my $id ( @$target_atom_id ) {
	@rotatable_path = [];
    }
}

1;
