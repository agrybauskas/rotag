package PDBxParser;

use strict;
use warnings;

use Exporter qw( import );
our @EXPORT_OK = qw( create_pdbx_entry
                     determine_residue_keys
                     extract
                     filter
                     filter_by_unique_residue_key
                     filter_new
                     identify_residue_atoms
                     indexed2raw
                     mark_selection
                     obtain_pdbx_data
                     obtain_pdbx_line
                     obtain_pdbx_loop
                     obtain_pdb_atom_site
                     obtain_atom_site
                     obtain_atom_sites
                     pdbx_indexed
                     pdbx_loop_to_array
                     pdbx_raw
                     raw2indexed
                     related_category_data
                     split_by
                     to_csv
                     to_pdbx
                     unique_residue_key );

use Carp;
use List::MoreUtils qw( any
                        uniq );
use Version qw( $VERSION );

our $VERSION = $VERSION;

# -------------------------------- Data structure ----------------------------- #

sub pdbx_raw
{
    my ( $pdbx_file, $data_identifier ) = @_;
    return obtain_pdbx_data( $pdbx_file, $data_identifier );
}

sub pdbx_indexed
{
    my ( $pdbx_file, $data_identifier, $options ) = @_;
    my $pdbx = obtain_pdbx_data( $pdbx_file, $data_identifier );
    raw2indexed( $pdbx, $options );
    return $pdbx;
}

# --------------------------------- PDBx parser ------------------------------- #

#
# Obtains pdbx data for the specified categories or items.
# Input:
#     $pdbx_file - PDBx file path;
#     $data_identifier - list of categories or items;
# Output:
#     %pdbx_data - data structure for pdbx data.
#

sub obtain_pdbx_data
{
    my ( $pdbx_file, $data_identifier, $options ) = @_;
    my ( $read_stream ) = ( $options->{'read_stream'} );

    $read_stream //= 0;

    my @pdbx_data = ();

    local @ARGV = ( $pdbx_file );

    if( defined $data_identifier && @{ $data_identifier } ) {
        my @pdbxs = ();

        # Slurp whole pdbx file.
        my $line_counter = 0;
        {
            local $/ = '';
            while( <> ) {
                push @pdbxs, map { 'data_' . $_ }
                             grep { $_ ne '' }
                             split /(?<!\S)data_/, $_;
                $line_counter++;
            }
        }

        warn "$pdbx_file is empty.\n" if $line_counter == 0;

        for my $pdbx ( @pdbxs ) {
            my %pdbx_data = ();
            %pdbx_data = (
                %pdbx_data,
                %{ obtain_pdbx_line( [ $pdbx ], $data_identifier ) } );
            %pdbx_data = (
                %pdbx_data,
                %{ obtain_pdbx_loop( [ $pdbx ], $data_identifier,
                                     { 'ignore_missing_categories' => 1 } ) } );

            push @pdbx_data, \%pdbx_data;
        }
    }

    if( $read_stream ) {
        return \@pdbx_data;
    } else {
        return $pdbx_data[0];
    }
}

#
# Obtains pdbx lines for a specified items.
# Input:
#     $pdbx - PDBx file content;
#     $items - list of desired items.
# Output:
#     %pdbx_line_data - hash of item values.
#

sub obtain_pdbx_line
{
    my ( $pdbx, $items ) = @_;

    my %pdbx_line_data;
    my %current_line_data;
    my $item_regexp = join q{|}, @{ $items };
    $item_regexp =~ s/\[/\\[/g;
    $item_regexp =~ s/\]/\\]/g;

    foreach( @{ $pdbx } ) {
        my %single_line_matches =
            ( m/($item_regexp|$item_regexp\.\S+)\s+(?!;)('.+'|\S+)/g );
        my %multi_line_matches =
            ( m/($item_regexp|$item_regexp\.\S+)\s+(\n;[^;]+;)/gx );
        %current_line_data = ( %single_line_matches, %multi_line_matches );
    }

    for my $key ( sort { $a cmp $b } keys %current_line_data ) {
        my ( $category, $attribute ) = split '\\.', $key;
        push @{$pdbx_line_data{$category}{'metadata'}{'attributes'}}, $attribute;
        push @{$pdbx_line_data{$category}{'data'}}, $current_line_data{$key};
        $pdbx_line_data{$category}{'metadata'}{'is_loop'} = 0;
        $pdbx_line_data{$category}{'metadata'}{'is_indexed'} = 0;
    }

    return \%pdbx_line_data;
}

#
# Obtains pdbx loops for a specified categories.
# Input:
#     $pdbx - PDBx file content;
#     $categories - list of specified categories.
# Output:
#     %pdbx_loop_data - data structure for loop data or list of data structure.
#

sub obtain_pdbx_loop
{
    my ( $pdbx, $categories, $options ) = @_;
    my ( $ignore_missing_categories ) = (
        $options->{'ignore_missing_categories'},
    );

    $ignore_missing_categories //= 0;

    my @categories;
    my @attributes;
    my @data; # Will be used for storing atom data temporarily.

    my $category_regexp = join q{|}, @{ $categories };
    $category_regexp =~ s/\[/\\[/g;
    $category_regexp =~ s/\]/\\]/g;
    my $is_reading_lines = 0; # Starts/stops reading lines at certain flags.

    my @pdbx = ( split /\n/, $pdbx->[0] );

    foreach( @pdbx ) {
        if( /^data_/ || ! @categories ) {
            push @categories, [];
            push @attributes, [];
            push @data, [];
        } elsif( /($category_regexp)[.](.+)$/x ) {
            if( !@{ $categories[-1] } || $categories[-1][-1] ne $1 ) {
                push @{ $categories[-1] }, $1;
                push @{ $attributes[-1] }, [];
                push @{ $data[-1] }, [];
            }
            push @{ $attributes[-1][-1] }, split q{ }, $2;
            $is_reading_lines = 1;
        } elsif( $is_reading_lines == 1 && /^_|loop_|#/ ) {
            if( $#categories eq $#{ $categories } ) { last; }
            $is_reading_lines = 0;
        } elsif( $is_reading_lines == 1 ) {
            my @current_data = ( $_ =~ m/('.+'|\S+)/g );
            push @{ $data[-1][-1] }, @current_data;
        }
    }

    # Checks the difference between the categories that were searched and
    # the ones that were found.
    if( ! $ignore_missing_categories ) {
        for my $current_categories ( @categories ) {
            for my $searched_category ( @{ $categories } ) {
                if( ! any { $searched_category eq $_ } @{ $current_categories }){
                    warn "'$searched_category' data was not found.\n";
                }
            }
        }
    }

    # Generates hash from three lists.
    my @pdbx_loop_data;
    for( my $i = 0; $i <= $#categories; $i++ ) {
        my %pdbx_loop_data;
        for( my $j = 0; $j <= $#{ $categories[-1] }; $j++ ) {
            $pdbx_loop_data{$categories[$i][$j]}{'metadata'}{'is_loop'} = 1;
            $pdbx_loop_data{$categories[$i][$j]}{'metadata'}{'is_indexed'} = 0;
            $pdbx_loop_data{$categories[$i][$j]}{'metadata'}{'attributes'} =
                $attributes[$i][$j];
            $pdbx_loop_data{$categories[$i][$j]}{'data'} =
                $data[$i][$j];
        }

        push @pdbx_loop_data, \%pdbx_loop_data;
    }

    return $pdbx_loop_data[0]; # FIXME: remove left-over data structure.
}

#
# Generates hash from pdbx loop data structure.
# Input:
#     $pdbx_loop_data - data structure (from obtain_pdbx_loop);
# Output:
#     @pdbx_loops - arrays of hashes data structure.
#

sub pdbx_loop_to_array
{
    my ( $pdbx_loop_data, $category ) = @_;

    my @attributes = @{ $pdbx_loop_data->{$category}{'metadata'}{'attributes'} };
    my @data = @{ $pdbx_loop_data->{$category}{'data'} };

    # Creates special data structure.
    my @pdbx_loops;
    my @data_row;
    my %data_row;

    my $attribute_count = scalar @attributes;
    my $data_count = scalar @data;

    for( my $pos = 0; $pos < $data_count; $pos += $attribute_count ) {
        @data_row = @{ data[$pos..$pos+$attribute_count-1] };
        %data_row = ();
        for( my $col = 0; $col <= $#data_row; $col++ ) {
            $data_row{$attributes[$col]} = $data_row[$col];
        }
        push @pdbx_loops, { %data_row };
    }

    return \@pdbx_loops;
}

#
# Creates a hash where for each unique residue key proper atom ids are assigned.
# Input:
#     $atom_site - $atom_site - atom data structure.
# Output:
#     %residue_atom_ids - hash of unique residue key and corresponding atom ids.
#

sub identify_residue_atoms
{
    my ( $atom_site, $options ) = @_;
    my ( $check_atom_names ) = $options->{'check_atom_names'};

    my $split_groups = split_by( { 'atom_site' => $atom_site,
                                   'attributes' => [ 'label_seq_id',
                                                     'label_asym_id',
                                                     'pdbx_PDB_model_num', ] } );

    my %residue_atom_ids;
    for my $atom_id ( keys %{ $atom_site } ) {
        my $unique_residue_key = unique_residue_key( $atom_site->{$atom_id} );
        my ( $residue_id, $residue_chain, $pdbx_model, $alt_id ) =
            split /,/, $unique_residue_key;
        my $split_group_entry = "${residue_id},${residue_chain},${pdbx_model}";
        my $related_residue_atom_ids = $split_groups->{$split_group_entry};

        # Splits related residues into alt id groups. Atom ids can be redundant.
        for my $related_atom_id ( @{ $related_residue_atom_ids } ) {
            my $related_alt_id =
                $atom_site->{$related_atom_id}{'label_alt_id'};

            if( $atom_id ne $related_atom_id ) {
                next if( $check_atom_names &&
                         $atom_site->{$atom_id}{'label_atom_id'} eq
                         $atom_site->{$related_atom_id}{'label_atom_id'} );

                if( $alt_id eq '.' ) {
                    push @{ $residue_atom_ids{$unique_residue_key} },
                        $related_atom_id;
                } elsif( $alt_id eq $related_alt_id || $related_alt_id eq '.') {
                    push @{ $residue_atom_ids{$unique_residue_key} },
                        $related_atom_id;
                }
            }
        }
    }

    return \%residue_atom_ids;
}

#
# Creates PDBx entry.
# Input:
#     $args - hash of all necessary attributes with corresponding values;
# Output:
#     PDBx STDOUT
#

sub create_pdbx_entry
{
    my ( $args ) = @_;
    my $atom_site = $args->{'atom_site'};
    my $atom_id = $args->{'id'};
    my $type_symbol = $args->{'type_symbol'};
    my $label_atom_id = $args->{'label_atom_id'};
    my $label_alt_id = $args->{'label_alt_id'};
    $label_alt_id //= q{.};
    my $label_comp_id = $args->{'label_comp_id'};
    my $label_asym_id = $args->{'label_asym_id'};
    my $label_entity_id = $args->{'label_entity_id'};
    $label_entity_id //= q{?};
    my $label_seq_id = $args->{'label_seq_id'};
    my $cartn_x = $args->{'cartn_x'};
    my $cartn_y = $args->{'cartn_y'};
    my $cartn_z = $args->{'cartn_z'};
    my $pdbx_model_num = $args->{'pdbx_PDB_model_num'};

    $atom_site->{$atom_id}{'group_PDB'} = 'ATOM';
    $atom_site->{$atom_id}{'id'} = $atom_id;
    $atom_site->{$atom_id}{'type_symbol'} = $type_symbol;
    $atom_site->{$atom_id}{'label_atom_id'} = $label_atom_id;
    $atom_site->{$atom_id}{'label_alt_id'} = $label_alt_id;
    $atom_site->{$atom_id}{'label_comp_id'} = $label_comp_id;
    $atom_site->{$atom_id}{'label_asym_id'} = $label_asym_id;
    $atom_site->{$atom_id}{'label_entity_id'} = $label_entity_id;
    $atom_site->{$atom_id}{'label_seq_id'} = $label_seq_id;
    $atom_site->{$atom_id}{'Cartn_x'} = $cartn_x;
    $atom_site->{$atom_id}{'Cartn_y'} = $cartn_y;
    $atom_site->{$atom_id}{'Cartn_z'} = $cartn_z;
    $atom_site->{$atom_id}{'pdbx_PDB_model_num'} = $pdbx_model_num;

    return;
}

sub related_category_data
{
    my ( $pdbx_data, $relationships ) = @_;

    my %related_category_data = ();

    for my $category ( sort keys %{ $relationships } ) {
        next if ! exists $pdbx_data->{$category};

        for my $related_category ( sort keys %{ $relationships->{$category} } ) {
            my @references = @{ $relationships->{$category}{$related_category} };
            for my $reference ( @references ) {
                my $keys = $reference->{'keys'};
                my $reference_keys = $reference->{'reference_keys'};
                $related_category_data{$category}{'keys'} = $reference_keys;
                $related_category_data{$category}{'reference_category'} =
                    $related_category;
                $related_category_data{$category}{'reference_keys'} =
                    $keys;

                raw2indexed( $pdbx_data,
                             { 'attributes' =>
                                   { $category => $reference_keys,
                                     $related_category => $keys },
                                   'is_unique' => 0 } );

                for my $key ( keys %{ $pdbx_data->{$related_category}{'data'} }){
                    $related_category_data{$category}{'data'}{$key} =
                        $pdbx_data->{$related_category}{'data'}{$key};
                }
            }
        }
    }

    return \%related_category_data;
}

# -------------------- Conversions between data structures -------------------- #

#
# Takes PDBx and indexes them by defined attributes.
# Input:
#     $pdbx - data structure generated with pdbx_raw().
#     $options->{attributes} - combination of attribute data that serves as
#     unique key.
#     $options->{'is_unique'} - checks the uniqueness.
# Output:
#     %indexed - returns indexed pdbx;
#

sub raw2indexed
{
    my ( $pdbx, $options ) = @_;
    my ( $attributes, $read_until_end, $default_is_unique ) = (
        $options->{'attributes'},
        $options->{'read_until_end'},
        $options->{'is_unique'},
    );

    $attributes //= {} unless $attributes;
    $default_is_unique //= 1;

    my %attributes = %{ $attributes };

    if( ! defined $pdbx || ! %{ $pdbx } ) {
        return {};
    }

    for my $category ( keys %{ $pdbx } ) {
        next if ! exists $pdbx->{$category} ||
                $pdbx->{$category}{'metadata'}{'is_indexed'} eq 1;

        my $keys = $attributes->{$category};
        my @attributes = @{ $pdbx->{$category}{'metadata'}{'attributes'} };
        my @data = @{ $pdbx->{$category}{'data'} };

        my $is_unique = $pdbx->{$category}{'metadata'}{'is_unique'};
        $is_unique //= $default_is_unique;

        # Creates special data structure.
        my @data_row;
        my %data_row;

        my $attribute_count = scalar @attributes;
        my $data_count = scalar @data;

        # Determines the positions of unique keys in attribute list.
        my @attribute_pos;
        for my $key ( @{ $keys } ) {
            for( my $i = 0; $i <= $#attributes; $i++ ) {
                if( $attributes[$i] eq $key ) {
                    push @attribute_pos, $i;
                    last;
                }
            }
        }

        my $auto_increment = 0;
        my %current_indexed = ();
        $current_indexed{'metadata'}{'attributes'} =
            $pdbx->{$category}{'metadata'}{'attributes'};
        $current_indexed{'metadata'}{'is_loop'} =
            $pdbx->{$category}{'metadata'}{'is_loop'};
        $current_indexed{'metadata'}{'is_indexed'} = 1;
        for( my $pos = 0; $pos <= $data_count - 1; $pos += $attribute_count ) {
            @data_row = @{ data[$pos..$pos+$attribute_count-1] };
            %data_row = ();
            for( my $col = 0; $col <= $#data_row; $col++ ) {
                $data_row{$attributes[$col]} = $data_row[$col];
            }

            my $current_key = join q{,}, map { $data_row[$_] } @attribute_pos;
            $current_key = $auto_increment unless $current_key;

            if( ! exists $current_indexed{'data'}{$current_key} ) {
                if( ! $is_unique ) {
                    $current_indexed{'data'}{$current_key} = [ { %data_row } ];
                } else {
                    $current_indexed{'data'}{$current_key} = { %data_row };
                }
            } else {
                if( ! $is_unique ) {
                    push @{ $current_indexed{'data'}{$current_key} },
                          { %data_row };
                } else {
                    confess 'unique key supposed to be unique.';
                }
            }

            $auto_increment++;
        }

        $pdbx->{$category} = { %current_indexed };
        $pdbx->{$category}{'metadata'}{'is_indexed'} = 1;
        $pdbx->{$category}{'metadata'}{'is_loop'} = 1;
    }

    return;
}

sub indexed2raw
{
    my ( $pdbx, $options ) = @_;
    my ( $categories, $attribute_order ) = (
        $options->{'categories'}, $options->{'attribute_order'}
    );

    $categories //= [ keys %{ $pdbx } ];

    for my $category ( @{ $categories } ) {
        next if ! exists $pdbx->{$category} ||
                $pdbx->{$category}{'metadata'}{'is_indexed'} eq 0;

        my $current_pdbx_indexed = $pdbx->{$category}{'data'};
        my $current_attribute_order = $attribute_order->{$category};

        $pdbx->{$category}{'metadata'}{'is_indexed'} = 0;
        $pdbx->{$category}{'data'} = (); # Resets the data.

        my @category_attributes =
            sort { $a cmp $b }
            uniq
            map { keys %{ $current_pdbx_indexed->{$_} } }
            keys %{ $current_pdbx_indexed };

        if( defined $current_attribute_order ) {
            @category_attributes = @{ $current_attribute_order };
        }

        $pdbx->{$category}{'metadata'}{'attributes'}=\@category_attributes;

        # HACK: should figure out how to deal with simple ids and combined
        # keys at the same time.
        for my $id ( sort { $a <=> $b } keys %{ $current_pdbx_indexed } ){
            for my $attribute ( @category_attributes ) {
                my $data_value = $current_pdbx_indexed->{$id}{$attribute};
                if( defined $data_value ) {
                    push @{ $pdbx->{$category}{'data'} }, $data_value;
                } else {
                    push @{ $pdbx->{$category}{'data'} }, '?';
                }
            }
        }
    }

    return;
}

# ----------------------------- Atom-site related ----------------------------- #

#
# From PDBx file, obtains data only from _atom_site category and outputs special
# data structure that represents atom data.
# Input:
#     $pdbx_file - PDBx file.
# Output:
#     %atom_site - indexed data structure.
#     E.g.: { 1 => { 'id' => 2,
#                    'label_atom_id' => 'CA',
#                    ... }
#             ... }
#

sub obtain_atom_site
{
    my ( $pdbx_file, $options ) = @_;
    my $atom_site =
        pdbx_indexed( $pdbx_file,
                      [ '_atom_site' ],
                      { 'attributes' => { '_atom_site' => [ 'id' ] } } );
    if( defined $atom_site && %{ $atom_site } ) {
        return $atom_site->{'_atom_site'}{'data'};
    } else {
        return {};
    }
}

#
# From PDBx file, obtains data only from _atom_site category and outputs special
# data structure that represents atom data. It is done for whole file or stream
# so, multiple data_ streams with _atom_site can be parsed.
# Input:
#     $pdbx_file - PDBx file.
# Output:
#     @atom_sites - list of special data structure.
#

sub obtain_atom_sites
{
    my ( $pdbx_file ) = @_;

    my $pdbxs = obtain_pdbx_data( $pdbx_file, [ '_atom_site' ],
                                  { 'read_stream' => 1 } );
    my @atom_sites = ();
    for my $pdbx ( @{ $pdbxs } ) {
        raw2indexed( $pdbx, { 'attributes' => { '_atom_site' => [ 'id' ] } } );
        push @atom_sites, $pdbx->{'_atom_site'}{'data'};
    }

    return \@atom_sites;
}

#
# Filters atom data structure according to specified attributes with include,
# exclude options.
# Input:
#     $args->{'atom_site'} - atom data structure;
#     $args->{'include'} - attribute selector that includes atom data structure.
#     Ex.:
#         { 'label_atom_id' => [ 'N', 'CA', 'CB', 'CD' ],
#           'label_comp_id' => [ 'A' ] };
#     $args->{'exclude'} - attribute selector that excludes atom data structure.
#     Selector data structure is the same as $args->{include};
#     $args->{'is_list'} - makes array instead of array of arrays;
#     $args->{'data_with_id'} - takes atom data structure and treats it as a
#     value and atom id - as a key;
#     $args->{'group_id'} - assigns the value of described group id.
# Output:
#     \%filtered_atoms- filtered atom data structure;
#

sub filter
{
    my ( $args ) = @_;
    my $atom_site = $args->{'atom_site'};
    my $include = $args->{'include'};
    my $exclude = $args->{'exclude'};
    my $data = $args->{'data'};
    my $is_list = $args->{'is_list'};
    my $data_with_id = $args->{'data_with_id'};
    my $group_id = $args->{'group_id'};

    if( ! defined $atom_site ) { confess 'no PDBx data structure was loaded '; }

    # Iterates through each atom in $atom_site and checks if atom specifiers
    # match up.
    my %filtered_atoms;

    # First, filters atoms that are described in $include specifier.
    if( defined $include && %{ $include } ) {
        for my $atom_id ( keys %{ $atom_site } ) {
            my $match_counter = 0; # Tracks if all matches occured.
            for my $attribute ( keys %{ $include } ) {
                if( exists $atom_site->{$atom_id}{$attribute} &&
                    any { $atom_site->{$atom_id}{$attribute} eq $_ }
                       @{ $include->{$attribute} } ) {
                    $match_counter += 1;
                } else {
                    last; # Terminates early if no match is found in specifier.
                }
            }
            if( $match_counter == scalar keys %{ $include } ) {
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
                if( exists $atom_site->{$atom_id}{$attribute} &&
                    any { $atom_site->{$atom_id}{$attribute} eq $_ }
                       @{ $exclude->{$attribute} } ) {
                    delete $filtered_atoms{$atom_id};
                    last;
                }
            }
        }
    }

    # TODO: again another iteration through atom data structure. Should look into
    # it how to reduce the quantity of iterations.
    if( defined $group_id ) {
        for my $atom_id ( keys %filtered_atoms ) {
            $filtered_atoms{$atom_id}{'[local]_selection_group'} = $group_id;
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
                    push @atom_data,
                         map { $filtered_atoms{$atom_id}{$_} } @{ $data };
                } else {
                    push @atom_data,
                         [ map { $filtered_atoms{$atom_id}{$_} } @{ $data } ];
                }
            }
            return \@atom_data;
        }
    }

    return \%filtered_atoms;
}

#
# Filters atom data structure according to specified attributes with include,
# exclude options.
# Input:
#     $options->{'include'} - attribute selector that includes atom data
#     structure.
#     Ex.:
#         { 'label_atom_id' => [ 'N', 'CA', 'CB', 'CD' ],
#           'label_comp_id' => [ 'A' ] };
#     $options->{'exclude'} - attribute selector that excludes atom data
#     $options->{'return_ids'} - flag that make function return filtered atom
#     ids instead of new AtomSite data structure.
#     $options->{'group_id'} - assigns the value of described group id to
#     '[local]_selection_group' attribute.
#     structure.
# Output:
#     \%filtered_atoms - filtered atom data structure;
#            or
#     \@filtered_atom_ids - list of filtered atom ids.
#

sub filter_new
{
    my ( $atom_site, $options ) = @_;
    my ( $include, $exclude, $return_data, $group_id, $selection_state ) =
        ( $options->{'include'}, $options->{'exclude'},
          $options->{'return_data'}, $options->{'group_id'},
          $options->{'selection_state'} );

    if( ! defined $atom_site ) {
        confess 'no atom were loaded to the AtomSite data structure';
    }

    $include //= {};
    $exclude //= {};

    # Generates hash maps for fast lookup.
    my %include_selector = ();
    for my $attribute ( keys %{ $include } ) {
        for my $value ( uniq @{ $include->{$attribute} } ) {
            $include_selector{$attribute}{$value} = 1;
        }
    }

    my %exclude_selector = ();
    for my $attribute ( keys %{ $exclude } ) {
        for my $value ( uniq @{ $exclude->{$attribute} } ) {
            $exclude_selector{$attribute}{$value} = 1;
        }
    }

    # Iterates through each atom in $self->{'atoms'} and checks if atom
    # specifiers match up.
    my %filtered_atoms;

    for my $atom_id ( keys %{ $atom_site } ) {
        my $keep_atom = 1;
        for my $attribute ( keys %{ $atom_site->{$atom_id} } ) {
            my $value = $atom_site->{$atom_id}{$attribute};
            if( exists $include_selector{$attribute} &&
                ( ! exists $include_selector{$attribute}{$value} ||
                  $include_selector{$attribute}{$value} != 1 ) ) {
                $keep_atom = 0;
                last;
            }

            if( exists $exclude_selector{$attribute} &&
                ( exists $exclude_selector{$attribute}{$value} &&
                  $exclude_selector{$attribute}{$value} == 1 ) ) {
                $keep_atom = 0;
                last;
            }
        }

        next if $keep_atom == 0;

        if( defined $group_id ) {
            $atom_site->{$atom_id}{'[local]_selection_group'} = $group_id;
        }

        if( defined $selection_state ) {
            $atom_site->{$atom_id}{'[local]_selection_state'} = $selection_state;
        }

        $filtered_atoms{$atom_id} = $atom_site->{$atom_id};
    }

    # Return object handle or atom ids depending on the flag.
    if( defined $return_data && $return_data eq 'id' ) {
        return [ keys %filtered_atoms ];
    } elsif( defined $return_data ) {
        return [ map { $atom_site->{$_}{$return_data} } keys %filtered_atoms ];
    } else {
        return \%filtered_atoms;
    }
}

#
# Filters atom site data structure by unique residue key.
# Input:
#     $atom_site - atom data structure.
#     $unique_residue_key - a composite key that identifies residue uniquely.
#     $include_dot - includes '.' alt id into the selection.
#     Ex.: '18,A,1,.'.
# Output:
#     %filtered_atoms - filtered atom data structure.
#

sub filter_by_unique_residue_key
{
    my ( $atom_site, $unique_residue_key, $include_dot ) = @_;
    my ( $residue_id, $residue_chain, $pdbx_model_num, $residue_alt ) =
        split /,/sxm, $unique_residue_key;
    my $filtered_atoms = filter( { 'atom_site' => $atom_site,
                                   'include' =>
                                   { 'label_seq_id' => [ $residue_id ],
                                     'label_asym_id' => [ $residue_chain ],
                                     'pdbx_PDB_model_num' => [ $pdbx_model_num ],
                                     'label_alt_id' =>
                                         [ $residue_alt,
                                           ( $include_dot ? '.' : () ) ] } } );
    return $filtered_atoms;
}

#
# Create unique residue key that consists of '_atom_site.label_seq_id',
# '_atom_site.label_asym_id', '_atom_site.pdbx_PDB_model_num' and
# '_atom_site.label_alt_id'.
# Input:
#     $atom - atom data structure.
# Output:
#     $unique_residue_key - unique residue key.
#

sub unique_residue_key
{
    my ( $atom ) = @_;
    return join q{,},
           map { $atom->{$_} }
               ( 'label_seq_id',
                 'label_asym_id',
                 'pdbx_PDB_model_num',
                 'label_alt_id', );
}

#
# TODO: maybe this function should be deprecated, because identify_residue_atoms
# does a better job.
# Generates a list of unique keys from the atom site.
# Input:
#     $atom_site - atom site data structure;
#     $options->{'exclude_dot'} - excludes label atom ids with '.' value, but
#     only if there are alternatives.
# Output:
#     @residue_unique_keys - list of determined unique keys.
#

sub determine_residue_keys
{
    my ( $atom_site, $options ) = @_;
    my ( $exclude_dot ) = $options->{'exclude_dot'};

    my @current_residue_unique_keys;
    for my $atom_id ( keys %{ $atom_site } ) {
        push @current_residue_unique_keys,
            unique_residue_key( $atom_site->{$atom_id} );
    }
    @current_residue_unique_keys = uniq @current_residue_unique_keys;

    my %residue_key_tree;
    for my $residue_unique_key ( @current_residue_unique_keys ) {
        my $reduced_unique_key = $residue_unique_key;
        my $alt_id = $residue_unique_key;
        $reduced_unique_key =~ s/^(.+,.+,.+),.+$/$1/g;
        $alt_id =~ s/^.+,.+,.+,(.+)$/$1/g;
        if( exists $residue_key_tree{$reduced_unique_key} ) {
            push @{ $residue_key_tree{$reduced_unique_key} },$residue_unique_key;
        } else {
            $residue_key_tree{$reduced_unique_key} = [ $residue_unique_key ];
        }
    }

    my @residue_unique_keys;
    for my $reduced_unique_key ( keys %residue_key_tree ) {
        my $residue_unique_keys = $residue_key_tree{$reduced_unique_key};
        if( scalar @{ $residue_unique_keys } > 1 ) {
            for my $i ( 0..$#{ $residue_unique_keys } ) {
                if( $exclude_dot && $residue_unique_keys->[$i] =~ m/\.$/ ) {
                    splice @{ $residue_unique_keys }, $i, 1;
                    last;
                }
            }
        }
        push @residue_unique_keys, @{ $residue_unique_keys };
    }

    return \@residue_unique_keys;
}

#
# Extracts information from atom data structure.
# Input:
#     $options{'data'} - list of data that should be extracted;
#     $options{'data_with_id'} - takes atom data structure and treats it as a
#     value and atom id - as a key;
#     $options{'is_list'} - - makes array instead of array of arrays;
# Output:
#

sub extract
{
    my ( $atom_site, $options ) = @_;
    my ( $data, $data_with_id, $is_list ) =
        ( $options->{'data'}, $options->{'data_with_id'},
          $options->{'is_list'} );

    my @atom_data;

    if( defined $data_with_id && $data_with_id ) {
        my %atom_data_with_id;

        # Simply iterates through $atom_site keys and extracts data
        # using data specifier and is asigned to atom id.
        for my $atom_id ( sort { $a <=> $b } keys %{ $atom_site } ) {
            $atom_data_with_id{$atom_id} =
                [ map { $atom_site->{$atom_id}{$_} } @{ $data } ];
        }

        return \%atom_data_with_id;
    } else {
        for my $atom_id ( sort { $a <=> $b } keys %{ $atom_site } ) {
            if( defined $is_list && $is_list ) {
                push @atom_data,
                    map { $atom_site->{$atom_id}{$_} } @{ $data };
            } else {
                push @atom_data,
                    [ map { $atom_site->{$atom_id}{$_} } @{ $data } ];
            }
        }

        return \@atom_data;
    }
}

#
# Splits up atom site to different groups by specified attributes.
# Input:
#     $args->{atom_site} - atom site data structure;
#     $args->{attributes} - list of attributes that atom site will be split by;
#     $args->{append_dot} - atoms that has 'label_alt_id' eq '.' to
#     corresponding groups.
# Output:
#     %split_groups - hash of atom site data structures.
#

sub split_by
{
    my ( $args ) = @_;
    my ( $atom_site, $attributes, $append_dot ) =
        ( $args->{'atom_site'}, $args->{'attributes'}, $args->{'append_dot'} );

    $attributes //=
        [ 'label_seq_id', 'label_asym_id', 'pdbx_PDB_model_num', 'label_alt_id'];
    $append_dot //= 0;

    my %split_groups;
    for my $atom_id ( sort keys %{ $atom_site } ) {
        # Creates group determining key that is used to sort atoms.
        my @attribute_values;
        for my $attribute ( @{ $attributes } ) {
            push @attribute_values, $atom_site->{$atom_id}{$attribute};
        }

        my $group_key = join q{,}, @attribute_values;

        if( exists $split_groups{$group_key} ) {
            push @{ $split_groups{$group_key} }, $atom_id;
        } else {
            $split_groups{$group_key} = [ $atom_id ];
        }
    }

    if( $append_dot ) {
        # Pre-determines position of attribute in unique key.
        my $alt_id_pos;
        for my $i ( 0..$#{ $attributes } ) {
            if( $attributes->[$i] eq 'label_alt_id' ) {
                $alt_id_pos = $i;
                last;
            }
        }

        # Defines relations alternative and origin atoms.
        my %unique_key_relations;
        my @origin_keys;
        for my $unique_key ( keys %split_groups ) {
            my @unique_key_attributes = split /,/sxm, $unique_key;
            if( $unique_key_attributes[$alt_id_pos] ne q{.} ) {
                $unique_key_attributes[$alt_id_pos] = '.';

                my $origin_key = join ',', @unique_key_attributes;
                push @origin_keys, $origin_key;

                $unique_key_relations{$unique_key} = $origin_key;
            }
        }

        # Appends origin atoms to alternative groups of atoms if necessary.
        for my $alt_key ( keys %unique_key_relations ) {
            if( exists $split_groups{$unique_key_relations{$alt_key}} ) {
                push @{ $split_groups{$alt_key} },
                     @{ $split_groups{$unique_key_relations{$alt_key}} };
            }
        }

        # Removes origin atoms that have alternatives.
        for my $origin_key ( uniq @origin_keys ) {
            delete $split_groups{$origin_key};
        }
    }

    return \%split_groups;
}

#
# Adds T (target), S (surrounding or selected), I (ignored) char to
# '[local]_selection_state' attribute data.
# Input:
#     $atom_site - atom site data structure;
#     $options{'target'} - list of ids of the target atom;
#     $options{'select'} - list of ids of the selected atom;
# Output:
#     adds markers to specified attribute field.
#

sub mark_selection
{
    my ( $atom_site, $options ) = @_;

    my ( $target_atom_ids, $selected_atom_ids ) =
        ( $options->{'target'}, $options->{'select'}, );

    for my $atom_id ( keys %{ $atom_site } ) {
        $atom_site->{$atom_id}{'[local]_selection_state'} = 'I';
    }

    for my $selected_atom_id ( @{ $selected_atom_ids } ) {
        $atom_site->{$selected_atom_id}{'[local]_selection_state'}='S';
    }

    for my $target_atom_id ( @{ $target_atom_ids } ) {
        $atom_site->{$target_atom_id}{'[local]_selection_state'} = 'T';
    }

    return;
}

# --------------------------------- PDB parser -------------------------------- #

sub obtain_pdb_atom_site
{
    my ( $pdbx_file ) = @_;

    my %atom_site = ();

    my $pdbx_PDB_model_num = 1;

    local @ARGV = ( $pdbx_file );
    while( <> ) {
        if( /^MODEL/ ) {
            my ( $model_id ) = m/^.{6}.{4}(.{4})$/;
            $pdbx_PDB_model_num = $model_id;
        } elsif( /^ATOM|^HETATM/ ) {
            # TODO: look if more '?!:' should be included in regex.
            my @atom_data = m/^(.{6})(.{5}).{1}(.{4})(.{1})(.{3}).{1}(.{1})(.{4}).{4}(.{8})(.{8})(.{8})(.{6})(.{6}).{10}(.{2})(?!:.{2})/;

            for my $atom_data_item ( @atom_data ) {
                $atom_data_item =~ s/ //g;
            }

            my ( $group_pdb, $id, $label_atom_id, $label_alt_id, $label_comp_id,
                 $label_asym_id, $label_seq_id, $cartn_x, $cartn_y, $cartn_z,
                 $occupancy, $B_iso_or_equiv_esd, $type_symbol,
                 $pdbx_formal_charge ) = @atom_data;

            $atom_site{$id}{'group_PDB'} = $group_pdb;
            $atom_site{$id}{'id'} = $id;
            $atom_site{$id}{'type_symbol'} = $type_symbol;
            $atom_site{$id}{'label_atom_id'} = $label_atom_id;
            $atom_site{$id}{'label_alt_id'} = '.' unless $label_alt_id;
            $atom_site{$id}{'label_comp_id'} = $label_comp_id;
            $atom_site{$id}{'label_asym_id'} = $label_asym_id;
            $atom_site{$id}{'label_seq_id'} = $label_seq_id;
            $atom_site{$id}{'Cartn_x'} = $cartn_x;
            $atom_site{$id}{'Cartn_y'} = $cartn_y;
            $atom_site{$id}{'Cartn_z'} = $cartn_z;
            $atom_site{$id}{'occupancy'} = $occupancy;
            $atom_site{$id}{'B_iso_or_equiv_esd'} = $B_iso_or_equiv_esd;
            $atom_site{$id}{'pdbx_formal_charge'} = $pdbx_formal_charge;
            $atom_site{$id}{'pdbx_PDB_model_num'} = $pdbx_PDB_model_num;
        }
    }

    return \%atom_site;
}

# --------------------------- Data structures to STDOUT ----------------------- #

#
# Converts pdbx data structure to PDBx.
# Input:
#     $options->{data_name} - data name of the PDBx;
#     $options->{category_order} - category list that should be included in the
#     output in certain order;
#     $options->{attribute_order} - attribute hash of lists that should be
#     included in the output in certain order;
#     $options->{add_attributes} - add list of attributes to existing data
#     structure;
#     $options->{tags} - tags that should be included in the output.
#     $options->{fh} - file handler.
# Output:
#     PDBx STDOUT.
#

sub to_pdbx
{
    my ( $pdbx_data, $options ) = @_;
    my ( $data_name, $categories, $category_order, $attribute_order, $tags,
         $add_attributes, $fh ) = (
        $options->{'data_name'},
        $options->{'categories'},
        $options->{'category_order'},
        $options->{'attribute_order'},
        $options->{'tags'},
        $options->{'add_attributes'},
        $options->{'fh'},
    );

    $data_name //= 'rotag';
    $category_order //= [
        '_atom_site', '_[local]_rotamer_angle', '_[local]_dihedral_angle',
        '_[local]_rotamer_energy', '_[local]_pairwise_energy', '_[local]_energy',
        '_[local]_rmsd'
    ];
    $attribute_order //= {
        '_atom_site' => [
            'group_PDB',
            'id',
            'type_symbol',
            'label_atom_id',
            'label_alt_id',
            'label_comp_id',
            'label_asym_id',
            'label_entity_id',
            'label_seq_id',
            'Cartn_x',
            'Cartn_y',
            'Cartn_z',
            'pdbx_PDB_model_num',
        ]
    };
    $fh //= \*STDOUT;

    my ( $current_categories ) = sort_by_list( [ sort keys %{ $pdbx_data } ],
                                               $category_order );
    $categories //= $current_categories;

    print {$fh} "data_${data_name}\n#\n";

    if( %{ $pdbx_data } ) {
        for my $category  ( @{ $categories } ) {
            next if defined $tags && ! any { $category eq $_ } @{ $tags };

            my $category_attribute_order = $attribute_order->{$category};
            if( ! defined $category_attribute_order ) {
                $category_attribute_order =
                    $pdbx_data->{$category}{'metadata'}{'attributes'};
            }

            my @append_attributes = ();
            for my $add_attribute ( @{ $add_attributes->{$category} } ) {
                if( ! any { $add_attributes eq $_ }
                         @{ $category_attribute_order } ) {
                    push @append_attributes, $add_attribute;
                }
            }
            push @{ $category_attribute_order }, @append_attributes;

            my ( undef, $current_attribute_order_table ) =
                sort_by_list( $category_attribute_order );

            if( $pdbx_data->{$category}{'metadata'}{'is_loop'} ) {
                print {$fh} "loop_\n";

                foreach( @{ $category_attribute_order } ) {
                    print {$fh} "$category.$_\n";
                }

                if( $pdbx_data->{$category}{'metadata'}{'is_indexed'} ) {
                    indexed2raw( $pdbx_data,
                                 { 'categories' => $categories,
                                   'attribute_order' => {
                                       $category => $category_attribute_order}});
                }

                my $attribute_array_length =
                    $#{ $pdbx_data->{$category}{'metadata'}{'attributes'} };
                my $data_array_length =
                    $#{ $pdbx_data->{$category}{'data'} };
                for( my $i = 0;
                     $i <= $data_array_length;
                     $i += $attribute_array_length + 1 ) {
                    my @current_data_list = ();
                    for my $j ( 0..$#{ $category_attribute_order } ) {
                        my $attribute = $category_attribute_order->[$j];
                        my $pos = $current_attribute_order_table->{$attribute};
                        if( defined $pos ) {
                            push @current_data_list,
                                $pdbx_data->{$category}{'data'}[$i+$pos];
                        } else {
                            push @current_data_list, '?'
                        }
                    }
                    print {$fh} join( q{ }, @current_data_list ), "\n" ;
                }
            } else { # PDBx line data.
                my @attributes =
                    @{ $pdbx_data->{$category}{'metadata'}{'attributes'} };
                my @data = @{ $pdbx_data->{$category}{'data'} };
                for( my $i = 0; $i <= $#attributes; $i++ ) {
                    printf {$fh} "%s.%s %s\n", $category, $attributes[$i],
                        $data[$i];
                }
            }

            print {$fh} "#\n";
        }
    }

    # return;
}

#
# Converts pdbx loop data structure to csv table.
# Input:
#     $pdbx_loops - data structure of pdbx_loops;
#     $attributes - columns that should be displayed.
# Output:
#     csv STDOUT
#

sub to_csv
{
    my ( $pdbx, $options ) = @_;
    my ( $category, $attributes ) = (
        $options->{'category'}, $options->{'attributes'}
    );

    $category //= (sort keys %{ $pdbx })[0]; # TODO: assigning through list
                                             # might loose the information.
    my $pdbx_raw = $pdbx->{$category};

    $attributes //= $pdbx_raw->{'metadata'}{'attributes'};

    if( defined $pdbx_raw ) {
        print {*STDOUT} join( ',', @{ $attributes } ), "\n";
    }

    my $attribute_array_length = $#{ $pdbx_raw->{'metadata'}{'attributes'} };
    my $data_array_length = $#{ $pdbx_raw->{'data'} };

    for( my $i = 0; $i <= $data_array_length; $i += $attribute_array_length + 1){
        print {*STDOUT} join( q{,}, @{ $pdbx_raw->{'data'} }
                                    [ $i..$i+$attribute_array_length ] ), "\n" ;
    }

    return;
}

sub sort_by_list
{
    my ( $unsorted_list, $sort_by_list ) = @_;

    $sort_by_list //= $unsorted_list;
    $unsorted_list //= $sort_by_list;

    my @sorted_list = sort @{ $unsorted_list };
    my %sort_order = ();

    my $upper_pos = $#{ $sort_by_list } + 1;
    for my $item ( @sorted_list ) {
        for my $pos ( 0..$#{ $sort_by_list } ) {
            if( $item eq $sort_by_list->[$pos] ) {
                $sort_order{$item} = $pos;
            }
        }

        # Move category to the back of the category list.
        if( ! exists $sort_order{$item} ) {
            $sort_order{$item} = $upper_pos;
            $upper_pos++;
        }
    }
    @sorted_list =
        sort { $sort_order{$a} <=> $sort_order{$b} } @sorted_list;

    return \@sorted_list, \%sort_order;
}

1;
