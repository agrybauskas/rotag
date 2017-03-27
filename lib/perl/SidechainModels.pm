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

{
my @visited_atoms;

sub follow_bonds
{
    my ( $atom_connections, $connection_template ) = @_;
    my %atom_connections = %$atom_connections;
    my %connection_template; # Used for tracking connections during recursion.

    if( ! $connection_template ) {
	%connection_template = %$atom_connections;
    }

    # print "....................\n";
    # print Dumper @visited_atoms;
    # print "-------------------\n";
    # print Dumper \%atom_connections;
    for my $atom ( keys %atom_connections ) {
    	if( ! grep { $atom eq $_ } @visited_atoms ) {
    	    push( @visited_atoms, $atom );

    	    if( scalar( @{ $atom_connections{$atom} } ) > 1 ) {

    		# Follows bond by creating new hash table from values of the
    		# current key.
    		my %current_connections;

    		foreach( @{ $connection_template{$atom} } ) {
    		    $current_connections{$_} = $connection_template{$_};
    		}

    		return follow_bonds( \%current_connections,
    				     \%connection_template );
    	    } elsif( scalar( @{ $atom_connections{$atom} } ) == 1 ) {
    		if( ! grep { $connection_template{$atom} } @visited_atoms ){
    		    push( @visited_atoms, $atom_connections{$atom}->[0] );

    		    my %current_connections;

    		    foreach( @{ $connection_template{$atom} } ) {
    			$current_connections{$_} = $connection_template{$_};
    		    }

    		    return follow_bonds( \%current_connections,
    					 \%connection_template );
    		}
    	    }
    	}
    }
    return @visited_atoms;
}
    # print "....................\n";
    # print Dumper @visited_atoms;
}

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
    follow_bonds( \%atom_connections );
}

1;
