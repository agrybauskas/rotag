package Utils;

use strict;
use warnings;

use Exporter qw( import );
our @EXPORT_OK = qw( select_atom_data );

use lib "../../lib/perl";
use PDBxParser qw( filter_atoms obtain_atom_site );

# -------------------------------- PDBx parsing ------------------------------- #

#
# Helps to shorten the code required to use PDBxParser subroutines.
#

#
# Returns atom data according to specified attribute values.
# Input:
#     $pdbx_file - PDBx file.
#     $atom_specifier - list of attributes and their values for specifying
#     desired atoms in string form.
#     Eg.: "label_atom_id CA,CB & label_comp_id SER"
#     $data_specifier - list of desired attribute values in string form.
#     Eg.: "Cartn_x,Cartn_y,Cartn_z".
# Output:
#     @atom_data - list of arrays containing specified attribute data.
#     Eg.: [ [ 3.0, 4.0, 2.0 ],
#            [ 3.4, 4.9, 2.5 ] ]
#

sub select_atom_data
{
    my ( $pdbx_file, $atom_specifier, $data_specifier,  ) = @_;

    my %atom_specifier = ( map { $_->[0] => [ split( ",", @$_[1] ) ] }
			   map { [ split( " ", $_ ) ] }
			   split( "&", $atom_specifier ) );
    my @data_specifier = split( ",", $data_specifier );

    my @atom_data =
    	@{ PDBxParser::select_atom_data(
    	       filter_atoms( obtain_atom_site( $pdbx_file ),
			     \%atom_specifier ),
    	       \@data_specifier ) };

    return \@atom_data;
}

1;
