package PseudoAtoms;

use strict;
use warnings;

use Exporter qw( import );
our @EXPORT_OK = qw( generate_pseudo );

use lib qw( ./ );
use CifParser qw( filter_atoms select_atom_data );
use SidechainModels qw( rotation_only );
use Data::Dumper;

# --------------------------- Generation of pseudo-atoms ---------------------- #

#
# Generates pseudo-atoms from side chain models that are written in equation
# form.
# Input  (3 arg): atom site data structure, hash of hashes for selecting atoms
#                 by attributes and hash of arrays that describe possible values
#                 of dihedral angles.
# Output (1 arg): atom site data structure with additional pseudo-atoms.
#
# Example of hash of arrays for describing dihedral angles.
# { "chi0" => [ [ 0, pi ], [ 1.5 * $pi, 2 * $pi ] ],
#   "chi1" => [ [ 0, 2 * $pi ] ] }
#

sub generate_pseudo
{
    my ( $atom_site, $atom_selector, $defined_angles, $resolution ) = @_;

    # Determines last id from set of atoms. It will be used for attaching id's
    # to generated pseudo-atoms.
    my $atom_ids = select_atom_data( [ "id" ], $atom_site );
    my @sorted_atom_ids = sort { $a <=> $b } map { $_->[0] } @{ $atom_ids };
    my $last_atom_id = $sorted_atom_ids[-1];

    # Generates model for selected atoms.
    $atom_site = rotation_only( $atom_site );

    my $target_atom_site = filter_atoms( $atom_selector, $atom_site );
    my $target_atom_ids = select_atom_data( [ "id" ], $target_atom_site );
    my @target_atom_ids = map { $_->[0] } @{ $target_atom_ids };

    for my $id ( @target_atom_ids ) {
	
    }
}

1;
