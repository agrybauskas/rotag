package CifParser;

use strict;
use warnings;

# ------------------------------ PDBx/mmCIF parser ---------------------------- #

#
# Extracts, filters and selects atom entries of PDBx/mmCIF files. "Attribute"
# coresponds to atom characteristics, such as atom or residues id, amino acid
# type and etc. Term "attribute" is used in CIF and mmCIF documentation.
#

#
# From mmCIF file, obtains data only from _atom_site category and outputs hash
# of compound data that represents attribute names and actual atom data.
# attribute names.
# Input  (1 arg): mmCIF file.
# Output (1 arg): hash of compound data structure.
#                 Ex.: { "attributes" => [ "group_id",
#                                          "id",
#                                          ...
#                                        ]
#                        "data" => { 1 => { "group_id" => "ATOM",
#                                           "id" => 1,
#                                           ...
#                                         }}}
#

sub obtain_atom_site
{
    my %atom_site;
    $atom_site{"attributes"} = [];
    $atom_site{"data"} = {};
    my @atom_data; # Will be used for storing of atom data temporarily.

    my $is_reading_lines = 0; # Starts/stops reading lines at certain flags.

    # Appends raw data to two keys in hash: attributes and data.
    foreach( @_ ) {
        if( $_ =~ /_atom_site\.(.+)\n$/ ) {
	    push( @{ $atom_site{"attributes"} }, split( " ", $1 ) );
            $is_reading_lines = 1;
        } elsif( $is_reading_lines == 1 && $_ =~ /^_|loop_|#/ ) {
            last;
        } elsif( $is_reading_lines == 1 ) {
            push( @atom_data, split( " ", $_ ) );
        }
    }

    # Converts atom_site data value from list to hash of hashes, that contain
    # attribute data assigned to actual values. ID attribute is used as key
    # accessing previously mentioned hashes.
    my @atom_data_row;
    my %atom_data_row;

    for( my $pos  = 0;
    	 $pos < $#atom_data;
    	 $pos += $#{ $atom_site{"attributes"} } + 1 ) {
	@atom_data_row =
	    @{ atom_data[$pos..$pos + $#{ $atom_site{"attributes"} }] };
	%atom_data_row = ();
	for( my $col = 0; $col <= $#atom_data_row; $col++ ) {
	    $atom_data_row{$atom_site{"attributes"}[$col]} =
		$atom_data_row[$col];
	}
	$atom_site{"data"}{$atom_data_row[1]} =
	    { %atom_data_row };
    }

    return \%atom_site;
}

#
# From mmCIF, extracts atoms with specified criteria, such as, atom type,
# residue id, chain id and etc.
# Input  (2 arg): array of hashes: atom specifier => value, mmCIF file converted
#                 to compound data using obtain_atom_site function.
# Output (1 arg): compound data, same as, obtain_atom_site output.
#

sub filter_atoms
{
    # Criteria for $atom_specifier:
    # ( "label_atom_id" => [ "SER" ],
    #   "label_comp_id" => [ "CA", "CB" ] ).
    my ( $atom_specifier, $atom_site ) = @_;

    my %filtered_atom_site;

    # Inherits list of attribute names from input.
    $filtered_atom_site{"attributes"} = $atom_site->{"attributes"};

    # Iterates through each atom in atom site and checks if atom specifiers
    # match up.
    my $match_counter; # Tracks if all matches occured.

    for my $atom ( keys %{ $atom_site->{"data"} } ) {
	$match_counter = 0;
	for my $attribute ( keys %$atom_specifier ) {
	    if( grep { $atom_site->{"data"}{$atom}{$attribute} eq $_ }
		@{ $atom_specifier->{$attribute} } ) {
		$match_counter += 1;
	    } else {
		last; # Terminates early if no match is found in any specifier.
	    }

	    if( $match_counter eq scalar( keys( %$atom_specifier ) ) ) {
		$filtered_atom_site{"data"}{$atom} = $atom_site->{"data"}{$atom};
	    }
	}
    }

    return \%filtered_atom_site;
}

#
# Returns specified attribute data in compound data form.
# Input  (2 arg): array of desired atom parameters, mmCIF file converted to
#                 compound data using obtain_atom_site function.
# Output (1 arg): compound data, same as, obtain_atom_site output.
#

sub select_atom_data
{
    # Extract only the data user wants.
    # E.g. of $data_specifier: ( "Cartn_x", "Cartn_y", "Cartn_z" ).
    my ( $data_specifier, $atom_site )  = @_;

    my @selected_atom_data;

    # Simply iterates through atom site keys and extracts data using data
    # specifier.
    for my $atom ( keys %{ $atom_site->{"data"} } ) {
	push( @selected_atom_data,
	      [ map { $atom_site->{"data"}{$atom}{$_} } @$data_specifier ] );
    }

    return \@selected_atom_data;
}

1;
