package PDBxParser;

use strict;
use warnings;

use Exporter qw( import );
our @EXPORT_OK = qw( filter_atoms
                     obtain_atom_site
                     select_atom_data
                     to_pdbx );
use Data::Dumper;
# --------------------------------- PDBx parser ------------------------------- #

#
# Extracts, filters and selects atom entries of PDBx files. "Attribute"
# coresponds to atom characteristics, such as atom or residue id, amino acid type
# and etc. Term "attribute" is used in PDBx documentation.
#

#
# From PDBx file, obtains data only from _atom_site category and outputs special
# data structure that represents atom data.
# Input:
#     $pdbx_file - PDBx file.
# Output:
#     %atom_site - special data structure.
#     Ex.: { 1 => { "group_id" => "ATOM",
#                   "id"       => 1,
#                   ... } }
#

sub obtain_atom_site
{
    my ( $pdbx_file ) = @_;

    my @atom_attributes;
    my @atom_data; # Will be used for storing atom data temporarily.
    my $is_reading_lines = 0; # Starts/stops reading lines at certain flags.

    # Finds atom data in PDBx and assigns to corresponding variables.
    open( my $fh, $pdbx_file );
    foreach( <$fh> ) {
        if( $_ =~ /_atom_site\.(.+)\n$/ ) {
    	    push( @atom_attributes, split( " ", $1 ) );
            $is_reading_lines = 1;
        } elsif( $is_reading_lines == 1 && $_ =~ /^_|loop_|#/ ) {
            last;
        } elsif( $is_reading_lines == 1 ) {
            push( @atom_data, split( " ", $_ ) );
        }
    }
    close( $fh );

    # Creates special data structure for describing atom site where atom id is
    # key in hash and hash value is hash describing atom data.
    my %atom_site;
    my @atom_data_row;
    my %atom_data_row;

    my $attribute_count = scalar( @atom_attributes );
    my $atom_data_count = scalar( @atom_data );

    for( my $pos = 0; $pos < $atom_data_count - 1; $pos += $attribute_count ) {
    	@atom_data_row =
    	    @{ atom_data[$pos..$pos+$attribute_count-1] };
    	%atom_data_row = ();
    	for( my $col = 0; $col <= $#atom_data_row; $col++ ) {
    	    $atom_data_row{$atom_attributes[$col]} =
    		$atom_data_row[$col];
    	}
    	$atom_site{$atom_data_row[1]} =
    	    { %atom_data_row };
    }

    return \%atom_site;
}

#
# From PDBx, extracts atoms with specified criteria, such as, atom type,
# residue id, chain id and etc.
# Input:
#     $atom_site - special data structure.
#     $atom_specifier - compound data structure for specifying desirable atoms.
#     Ex. { "label_atom_id" => [ "SER" ],
#           "label_comp_id" => [ "CA", "CB" ] }
# Output:
#     %filtered_atoms - filtered special data structure.
#

sub filter_atoms
{
    my ( $atom_site, $atom_specifier ) = @_;

    # Produces attribute list from subroutine input - $atom_site, by catching all
    # mentioned attributes.
    my @atom_attributes;

    for my $atom_id ( keys %{ $atom_site } ) {
	for my $attribute ( keys $atom_site{$atom_id} ) {
	    foreach ( @atom_attributes ) { # Checks, if attribute was already
	    }
	}
    }

    # $filtered_atom_site{"attributes"} = $atom_site->{"attributes"};
    # my %filtered_atoms;
    # # Iterates through each atom in atom site and checks if atom specifiers
    # # match up.
    # my $match_counter; # Tracks if all matches occured.

    # for my $atom ( keys %{ $atom_site->{"data"} } ) {
    # 	$match_counter = 0;
    # 	for my $attribute ( keys %$atom_specifier ) {
    # 	    if( exists $atom_site->{"data"}{$atom}{$attribute}
    # 	     && grep { $atom_site->{"data"}{$atom}{$attribute} eq $_ }
    # 		@{ $atom_specifier->{$attribute} } ) {
    # 		$match_counter += 1;
    # 	    } else {
    # 		last; # Terminates early if no match is found in any specifier.
    # 	    }
    # 	    if( $match_counter eq scalar( keys( %$atom_specifier ) ) ) {
    # 		$filtered_atom_site{"data"}{$atom} = $atom_site->{"data"}{$atom};
    # 	    }
    # 	}
    # }

    # return \%filtered_atom_site;
}

#
# Returns specified attribute data in compound data form.
# Input  (2 arg): array of desired atom parameters, mmCIF/PDBx file converted to
#                 compound data using obtain_atom_site function.
#                 E.g. of $data_specifier:
#                    ( "Cartn_x", "Cartn_y", "Cartn_z" ).
# Output (1 arg): compound data, same as, obtain_atom_site output.
#

sub select_atom_data
{
    # Extract only the data user wants.
    my ( $data_specifier, $atom_site )  = @_;

    my @selected_atom_data;

    # Simply iterates through atom site keys and extracts data using data
    # specifier.
    for my $atom ( sort { $a <=> $b } keys( %{ $atom_site->{"data"} } ) ) {
	push( @selected_atom_data,
	      [ map { $atom_site->{"data"}{$atom}{$_} } @{ $data_specifier } ] );
    }

    return \@selected_atom_data;
}

# --------------------------- Data structure to STDOUT ------------------------ #

#
# Converts special atom site data structure back to PDBx, XYZ (for Jmol) and etc.
#

#
# Converts to truncated PDBx format file.
# Input:
#     $atom_site - special data structure.
#     $atom_attributes - list of attributes
# Output:
#     STDOUT - PDBx file.

sub to_pdbx
{
    my ( $atom_site, $data_name, $atom_attributes ) = @_;

    # Assigns default data name and attributes if $atom_attributes variables are
    # undefined.
    $data_name //= "testing";
    $atom_attributes //= [ "group_PDB",
			   "id",
			   "type_symbol",
			   "label_atom_id",
			   "label_alt_id",
			   "label_comp_id",
			   "label_asym_id",
			   "label_entity_id",
			   "label_seq_id",
			   "pdbx_PDB_ins_code",
			   "Cartn_x",
			   "Cartn_y",
			   "Cartn_z",
			   "occupancy",
			   "B_iso_or_equiv",
			   "Cartn_x_esd",
			   "Cartn_y_esd",
			   "Cartn_z_esd",
			   "occupancy_esd",
			   "B_iso_or_equiv_esd",
			   "pdbx_formal_charge",
			   "auth_seq_id",
			   "auth_comp_id",
			   "auth_asym_id",
			   "auth_atom_id",
			   "pdbx_PDB_model_num" ];

    # Sends PDBx to STDOUT.
    print "data_$data_name\n";
    print "loop_\n";

    for my $attribute ( @{ $atom_attributes } ) {
    	print "_atom_site.$attribute\n";
    }

    for my $id ( sort { $a <=> $b } keys %{ $atom_site } ) {
    	for( my $i = 0; $i <= $#{ $atom_attributes }; $i++ ) {
    	    if( $i % ( $#{ $atom_attributes } + 1) != 0 ) {
    		print( " ", $atom_site->{$id}{$atom_attributes->[$i]}, " " );
    	    } else {
    		print( "\n", $atom_site->{$id}{$atom_attributes->[$i]} );
    	    }
    	}
    }
    print( "\n" );
}

1;
