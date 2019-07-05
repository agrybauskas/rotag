package PDBxParser_new;

use strict;
use warnings;

use Exporter qw( import );
our @EXPORT_OK = qw( indexed2raw
                     obtain_pdbx_data
                     pdbx_indexed
                     pdbx_raw
                     raw2indexed
                     to_pdbx
                     related_category_data );

use Carp;
use List::MoreUtils qw( any
                        uniq );
use Version qw( $VERSION );

our $VERSION = $VERSION;

# --------------------------------- PDBx parser ------------------------------- #

sub pdbx_raw
{
    my ( $pdbx_file, $data_identifier ) = @_;
    return obtain_pdbx_data( $pdbx_file, $data_identifier );
}

sub pdbx_indexed
{
    my ( $pdbx_file, $data_identifier, $options ) = @_;
    return raw2indexed( obtain_pdbx_data( $pdbx_file, $data_identifier ),
                        $options );
}

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
    my ( $pdbx_file, $data_identifier ) = @_;
    my %pdbx_data = ();

    local @ARGV = ( $pdbx_file );

    if( defined $data_identifier && @{ $data_identifier } ) {
        my @pdbx = ();

        # Slurp whole pdbx file.
        {
            local $/ = '';
            while( <> ) {
                push @pdbx, $_;
            }
        }

        %pdbx_data = (
            %pdbx_data,
            %{ obtain_pdbx_line( \@pdbx, $data_identifier ) } );
        %pdbx_data = (
            %pdbx_data,
            %{ obtain_pdbx_loop( \@pdbx, $data_identifier,
                                 { 'ignore_missing_categories' => 1 } ) } );
    }

    return \%pdbx_data;
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
#     $options->{'read_until_end'} - reads whole pdbx file or stdin.
# Output:
#     %pdbx_loop_data - data structure for loop data or list of data structure.
#

sub obtain_pdbx_loop
{
    my ( $pdbx, $categories, $options ) = @_;
    my ( $read_until_end, $ignore_missing_categories ) = (
        $options->{'read_until_end'},
        $options->{'ignore_missing_categories'},
    );

    $read_until_end //= 0;
    $ignore_missing_categories //= 0;

    my @categories;
    my @attributes;
    my @data; # Will be used for storing atom data temporarily.

    my $category_regexp = join q{|}, @{ $categories };
    $category_regexp =~ s/\[/\\[/g;
    $category_regexp =~ s/\]/\\]/g;
    my $is_reading_lines = 0; # Starts/stops reading lines at certain flags.

    my @pdbx = ( split /\n/, $pdbx->[0] );
    my $line_counter = 0;

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
            if( $#categories eq $#{ $categories } && ! $read_until_end ) { last; }
            $is_reading_lines = 0;
        } elsif( $is_reading_lines == 1 ) {
            my @current_data = ( $_ =~ m/('.+'|\S+)/g );
            push @{ $data[-1][-1] }, @current_data;
        }
        $line_counter++;
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

    warn "pdbx file is empty.\n"
        if $line_counter == 0 && ! $ignore_missing_categories;

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

    return $pdbx_loop_data[0] if ! $read_until_end;
    return \@pdbx_loop_data;
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

                my $category_data =
                    raw2indexed( { $related_category =>
                                       $pdbx_data->{$related_category} },
                                 { 'attributes' =>
                                       { $related_category => $keys },
                                   'is_unique' => 0 } );
                my $related_category_data =
                    raw2indexed( { $category => $pdbx_data->{$category} },
                                 { 'attributes' =>
                                       { $category => $reference_keys },
                                   'is_unique' => 0 } );
                for my $key ( keys %{ $category_data->{$related_category}{'data'} } ) {
                    $related_category_data{$category}{'data'}{$key} =
                        $related_category_data->{$category}{'data'}{$key};
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
#     $attributes - combination of attribute data that serves as unique key.
#     $options->{'read_until_end'} - reads whole pdbx file or stdin.
# Output:
#     %indexed - returns indexed pdbx;
#

sub raw2indexed
{
    my ( $pdbx_raw, $options ) = @_;
    my ( $attributes, $read_until_end, $default_is_unique ) = (
        $options->{'attributes'},
        $options->{'read_until_end'},
        $options->{'is_unique'},
    );

    $attributes //= {} unless $attributes;
    $read_until_end //= 0;
    $default_is_unique //= 1;

    my %attributes = %{ $attributes };

    if( ! defined $pdbx_raw || ! %{ $pdbx_raw } ) {
        return {};
    }

    my %indexed = ();

    for my $category ( keys %{ $pdbx_raw } ) {
        my $keys = $attributes->{$category};
        my @attributes = @{ $pdbx_raw->{$category}{'metadata'}{'attributes'} };
        my @data = @{ $pdbx_raw->{$category}{'data'} };

        my $is_unique = $pdbx_raw->{$category}{'metadata'}{'is_unique'};
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
            $pdbx_raw->{$category}{'metadata'}{'attributes'};
        $current_indexed{'metadata'}{'is_loop'} =
            $pdbx_raw->{$category}{'metadata'}{'is_loop'};
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

        $indexed{$category} = { %current_indexed };
    }

    return \%indexed;
}

sub indexed2raw
{
    my ( $pdbx_indexed ) = @_;

    my %pdbx_raw = ();
    for my $category ( keys %{ $pdbx_indexed } ) {
        my $current_pdbx_indexed = $pdbx_indexed->{$category}{'data'};

        $pdbx_raw{$category}{'metadata'}{'is_loop'} =
            $pdbx_indexed->{$category}{'metadata'}{'is_loop'};
        $pdbx_raw{$category}{'metadata'}{'is_indexed'} = 0;


        my ( $first_pdbx_id ) = sort keys %{ $current_pdbx_indexed };
        if( ref $current_pdbx_indexed->{$first_pdbx_id} eq 'HASH' ) {
            my @category_attributes =
                sort { $a cmp $b }
                uniq
                map { keys %{ $current_pdbx_indexed->{$_} } }
                keys %{ $current_pdbx_indexed };
            $pdbx_raw{$category}{'metadata'}{'attributes'}=\@category_attributes;

            # HACK: should figure out how to deal with simple ids and combined
            # keys at the same time.
            for my $id ( sort { $a <=> $b } keys %{ $current_pdbx_indexed } ){
                for( my $i = 0; $i <= $#category_attributes; $i++ ) {
                    my $data_value =
                        $current_pdbx_indexed->{$id}{$category_attributes[$i]};
                    if( defined $data_value ) {
                        push @{ $pdbx_raw{$category}{'data'} }, $data_value;
                    } else {
                        push @{ $pdbx_raw{$category}{'data'} }, '?';
                    }
                }
            }
        } elsif( ref $current_pdbx_indexed->{$first_pdbx_id} eq 'ARRAY' ) {
            my @category_attributes =
                sort { $a cmp $b }
                uniq
                map { keys %{ $_ } }
                map { @{ $current_pdbx_indexed->{$_} } }
                keys %{ $current_pdbx_indexed };
            $pdbx_raw{$category}{'metadata'}{'attributes'}=\@category_attributes;

            for my $id ( sort { $a cmp $b } keys %{ $current_pdbx_indexed } ){
                for my $group ( @{ $current_pdbx_indexed->{$id} } ) {
                    for( my $i = 0; $i <= $#category_attributes; $i++ ) {
                        my $data_value = $group->{$category_attributes[$i]};
                        if( defined $data_value ) {
                            push @{ $pdbx_raw{$category}{'data'} }, $data_value;
                        } else {
                            push @{ $pdbx_raw{$category}{'data'} }, '?';
                        }
                    }
                }
            }
        }
    }

    return \%pdbx_raw;
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
    my ( $data_name, $category_order, $attribute_order, $tags, $add_attributes,
         $fh ) = (
        $options->{'data_name'},
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
    $add_attributes //= { '_atom_site' => [ '_[local]_selection_state',
                                            '_[local]_selection_group' ] };
    $fh //= \*STDOUT;

    my ( $categories ) =
        sort_by_list( [ sort keys %{ $pdbx_data } ], $category_order );

    print {$fh} "data_${data_name}\n#\n";

    if( %{ $pdbx_data } ) {
        for my $category  ( @{ $categories } ) {
            next if defined $tags &&
                  ! any { $category eq $_ } @{ $tags };

            my $category_attribute_order = $attribute_order->{$category};
            if( ! defined $category_attribute_order ) {
                $category_attribute_order = $pdbx_data->{$category}{'metadata'}
                                                                   {'attributes'};
            }

            my @append_attributes = ();
            for my $add_attribute ( @{ $add_attributes->{$category} } ) {
                if( ! any { $add_attributes eq $_ }
                         @{ $category_attribute_order } ) {
                    push @append_attributes, $add_attribute;
                }
            }
            push @{ $category_attribute_order }, @append_attributes;

            my ( undef, $current_attribute_order ) =
                sort_by_list( $category_attribute_order,
                              $pdbx_data->{$category}{'metadata'}{'attributes'});

            if( $pdbx_data->{$category}{'metadata'}{'is_loop'} ) {
                print {$fh} "loop_\n";

                foreach( @{ $category_attribute_order } ) {
                    print {$fh} "$category.$_\n";
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
                        my $pos = $current_attribute_order->{$attribute};
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

    return;
}

sub sort_by_list
{
    my ( $unsorted_list, $sort_by_list ) = @_;

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
