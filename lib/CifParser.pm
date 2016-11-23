package CifParser;

use strict;
use warnings;

use List::MoreUtils qw( first_index );     # TODO: remove and replace dependency.

no if $] >= 5.017011,                       # WARNING: in newer versions of Perl,
    warnings => 'experimental::smartmatch'; # smartmach became experimental.

# ------------------------------ PDBx/mmCIF parser ---------------------------- #

#
# In this block of code, functions extract, filter and select atom entries of
# PDBx/mmCIF files. "Attribute" coresponds to atom characteristics, such as
# atom or residues id, amino acid type and etc. Term "attribute" is used in CIF
# and mmCIF documentation.
#

#
# From mmCIF file, obtains data only from _atom_site category and outputs 1x2
# array of attribute names and attribute data respectively.
#

sub obtain_atom_site
{
    my @atom_attributes;
    my @atom_data;

    my $is_reading_lines = 0;      # Starts/stops reading lines at certain flags.

    foreach( @_ ) {
        if( $_ =~ /_atom_site\.(.+)\n$/ ) {
            push( @atom_attributes, split( " ", $1 ) );
            $is_reading_lines = 1;
        } elsif( $is_reading_lines == 1 && $_ =~ /^_|loop_|#/ ) {
            last;
        } elsif( $is_reading_lines == 1 ) {
            push( @atom_data, split( " ", $_ ) );
        }
    }

    return \@atom_attributes, \@atom_data;
}

#
# From mmCIF file, extracts atoms with specified criteria, such as, atom type,
# residue id, chain id and etc.
#

sub filter_atoms
{
    # Criteria for desirable atoms. 
    # E.g. [ "label_atom_id" => ["SER"], 
    #        "label_atom_id" => ["CA", "CB"] ].
    my $atom_specifiers = shift;
    my %atom_specifiers = @$atom_specifiers;
    my @mmcif_stdin     = @_;

    my @atom_site = obtain_atom_site( @mmcif_stdin );
    my @atom_attributes = @{ $atom_site[0] };
    my @atom_data = @{ $atom_site[1] };

    my @attribute_pos; # The position of specified atom attributes in actual
                       # list of attributes of CIF file.

    for my $attribute ( keys %atom_specifiers ) {
        if( $attribute ~~ @atom_attributes ) {
            push( @attribute_pos,
                  first_index{ $_ eq $attribute } @atom_attributes );
        }
    }

    my @filtered_atoms;

    my @atom_data_row;
    my @spec_attributes;
    my @specified_data;

    for( my $pos  = 0; $pos < $#atom_data; $pos += $#atom_attributes + 1) {
        @atom_data_row = @{ atom_data[$pos..$pos + $#atom_attributes] };
        @spec_attributes = map { $atom_data_row[$_] } @attribute_pos;
        @specified_data  = map { $atom_specifiers{$_} }
                           map { $atom_attributes[$_] } @attribute_pos;

        if( @spec_attributes ~~ @specified_data ) {
            push( @filtered_atoms, @atom_data_row );
        }
    }

    return \@atom_attributes, \@filtered_atoms;
}

#
# Returns specified attribute data.
#

sub select_atom_data
{
    my $atom_specifiers = shift;
    my %atom_specifiers = @$atom_specifiers;
    my @data_specifier  = shift;   # Extract only the data user wants. 
                                   # E.g. [ "Cartn_x", "Cartn_y", "Cartn_z" ].
    my @mmcif_stdin = @_; # mmCIF type data.

    my @filtered_atoms = filter_atoms( $atom_specifiers, @mmcif_stdin );
    my @attribute_data = @{ $filtered_atoms[0] };
    my @atom_data = @{ $filtered_atoms[1] };

    my @attribute_pos;

    for my $attribute ( @{ $data_specifier[0] } ) {
        if( $attribute ~~ @attribute_data ) {
            push( @attribute_pos,
                  first_index{ $_ eq $attribute } @attribute_data );
        }
    }

    my @atom_data_row;
    my @selected_atom_data;

    for( my $pos  = 0; $pos < $#atom_data; $pos += $#attribute_data + 1) {
	@atom_data_row = @{ atom_data[$pos..$pos + $#attribute_data] };
	push( @selected_atom_data,
	      [ map { $atom_data_row[$_] } @attribute_pos ] );
    }

    return \@selected_atom_data;
}

1;
