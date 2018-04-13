package PDBxParser;

use strict;
use warnings;

use Exporter qw( import );
our @EXPORT_OK = qw( create_pdbx_entry
                     filter
                     obtain_atom_site
                     obtain_category_data
                     to_pdbx );

# --------------------------------- PDBx parser ------------------------------- #

#
# Extracts, filters and selects atom entries of PDBx files. "Attribute"
# coresponds to atom characteristics, such as atom or residue id, amino acid type
# and etc. Term "attribute" is used in PDBx documentation.
#

sub obtain_category_data
{
    my ( $pdbx_file, $categories ) = @_;

    my @category_groups;
    my @categories;
    my @attributes;

    my $is_reading_lines = 0; # Starts/stops reading lines at certain flags.
    my $is_loop = 0;
    my $is_multiline = 0;

    my $regexp_pattern = join( "|", @{ $categories } );

    @ARGV = ( $pdbx_file );
    while( <> ) {
	if( $_ =~ /($regexp_pattern)\.(.+)\n$/x ) {

	}
	#     push( @atom_attributes, split( " ", $1 ) );
	#     $is_reading_lines = 1;
	# } elsif( $is_reading_lines == 1 && $_ =~ /^_|loop_|#/ ) {
	#     last;
	# } elsif( $is_reading_lines == 1 ) {
	#     push( @atom_data, split( " ", $_ ) );
	# }
    }
}

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

    @ARGV = ( $pdbx_file );
    while( <> ) {
	if( $_ =~ /_atom_site\.(.+)\n$/x ) {
	    push( @atom_attributes, split( " ", $1 ) );
	    $is_reading_lines = 1;
	} elsif( $is_reading_lines == 1 && $_ =~ /^_|loop_|#/ ) {
	    last;
	} elsif( $is_reading_lines == 1 ) {
	    push( @atom_data, split( " ", $_ ) );
	}
    }

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

sub filter
{
    my ( $args ) = @_;
    my $atom_site = $args->{"atom_site"};
    my $include = $args->{"include"};
    my $exclude = $args->{"exclude"};
    my $data = $args->{"data"};
    my $is_list = $args->{"is_list"};
    my $data_with_id = $args->{"data_with_id"};

    die( "No PDBx data structure was loaded " ) if ! defined $atom_site;

    # Iterates through each atom in $atom_site and checks if atom specifiers
    # match up.
    my %filtered_atoms;

    # First, filters atoms that are described in $include specifier.
    if( defined $include && %{ $include } ) {
    	for my $atom_id ( keys %{ $atom_site } ) {
    	    my $match_counter = 0; # Tracks if all matches occured.
    	    for my $attribute ( keys %{ $include } ) {
    		if( exists $atom_site->{$atom_id}{$attribute}
    	         && grep { $atom_site->{$atom_id}{$attribute} eq $_ }
    		    @{ $include->{$attribute} } ) {
    		    $match_counter += 1;
    		} else {
    		    last; # Terminates early if no match is found in specifier.
    		}
    	    }
	    if( $match_counter == scalar( keys %{ $include } ) ) {
	    	$filtered_atoms{$atom_id} = $atom_site->{$atom_id};
	    }
    	}
    } else {
    	%filtered_atoms = %{ $atom_site };
    }

    # Then filters out atoms that are in $exclude specifier.
    if( defined $exclude && %{ $exclude } ) {
    	for my $atom_id ( keys %filtered_atoms ) {
    	    for my $attribute ( keys %{ $exclude } ) {
    	    	if( exists $atom_site->{$atom_id}{$attribute}
    	         && grep { $atom_site->{$atom_id}{$attribute} eq $_ }
    	    	    @{ $exclude->{$attribute} } ) {
    	    	    delete $filtered_atoms{$atom_id};
    		    last;
    	    	}
    	    }
    	}
    }

    # Extracts specific data, if defined.
    if( defined $data && @{ $data } ) {
    	# Simply iterates through $atom_site keys and extracts data using data
    	# specifier.
    	my @atom_data;
    	if( defined $data_with_id && $data_with_id ) {
    	    my %atom_data_with_id;
    	    # Simply iterates through $atom_site keys and extracts data using
    	    # data specifier and is asigned to atom id.
    	    for my $atom_id ( sort { $a <=> $b } keys %{ $atom_site } ) {
    		$atom_data_with_id{$atom_id} =
    		    [ map { $atom_site->{$atom_id}{$_} } @{ $data } ];
    	    }
    	    return \%atom_data_with_id;

    	} else {
    	    for my $atom_id ( sort { $a <=> $b } keys %filtered_atoms ) {
    		if( defined $is_list && $is_list ) {
    		    push( @atom_data,
    			  map { $filtered_atoms{$atom_id}{$_} } @{ $data } );
    		} else {
    		    push( @atom_data,
    			  [ map { $filtered_atoms{$atom_id}{$_} } @{ $data } ] );
    		}
    	    }
    	    return \@atom_data;
    	}
    }

    return \%filtered_atoms;
}

sub create_pdbx_entry
{
    my ( $args ) = @_;
    my $atom_site = $args->{"atom_site"};
    my $atom_id = $args->{"id"};
    my $type_symbol = $args->{"type_symbol"};
    my $label_atom_id = $args->{"label_atom_id"};
    my $label_alt_id = $args->{"label_alt_id"};
    $label_alt_id //= ".";
    my $label_comp_id = $args->{"label_comp_id"};
    my $label_asym_id = $args->{"label_asym_id"};
    my $label_entity_id = $args->{"label_entity_id"};
    $label_entity_id //= "?";
    my $label_seq_id = $args->{"label_seq_id"};
    my $cartn_x = $args->{"cartn_x"};
    my $cartn_y = $args->{"cartn_y"};
    my $cartn_z = $args->{"cartn_z"};
    my $auth_seq_id = $args->{"auth_seq_id"};
    my $auth_comp_id = $label_comp_id;
    my $auth_asym_id = $label_asym_id;
    my $auth_atom_id = $label_atom_id;
    my $pdbx_PDB_model_num = $label_entity_id;

    $atom_site->{$atom_id}{"group_PDB"} = "ATOM";
    $atom_site->{$atom_id}{"id"} = $atom_id;
    $atom_site->{$atom_id}{"type_symbol"} = $type_symbol;
    $atom_site->{$atom_id}{"label_atom_id"} = $label_atom_id;
    $atom_site->{$atom_id}{"label_alt_id"} = $label_alt_id;
    $atom_site->{$atom_id}{"label_comp_id"} = $label_comp_id;
    $atom_site->{$atom_id}{"label_asym_id"} = $label_asym_id;
    $atom_site->{$atom_id}{"label_entity_id"} = $label_entity_id;
    $atom_site->{$atom_id}{"label_seq_id"} = $label_seq_id;
    $atom_site->{$atom_id}{"Cartn_x"} = $cartn_x;
    $atom_site->{$atom_id}{"Cartn_y"} = $cartn_y;
    $atom_site->{$atom_id}{"Cartn_z"} = $cartn_z;
    $atom_site->{$atom_id}{"auth_seq_id"} = $auth_seq_id;
    $atom_site->{$atom_id}{"auth_comp_id"} = $auth_comp_id;
    $atom_site->{$atom_id}{"auth_asym_id"} = $auth_asym_id;
    $atom_site->{$atom_id}{"auth_atom_id"} = $auth_atom_id;
    $atom_site->{$atom_id}{"pdbx_PDB_model_num"} = $pdbx_PDB_model_num;

    return;
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
		if( exists $atom_site->{$id}{$atom_attributes->[$i]} ) {
		    print( " ", $atom_site->{$id}{$atom_attributes->[$i]}, " " );
		} else {
		    print( " ? " );
		}
    	    } else {
		if( exists $atom_site->{$id}{$atom_attributes->[$i]} ) {
		    print( "\n", $atom_site->{$id}{$atom_attributes->[$i]} );
		} else {
		    print( "\n." );
		}
    	    }
    	}
    }
    print( "\n" );

    return;
}

1;
