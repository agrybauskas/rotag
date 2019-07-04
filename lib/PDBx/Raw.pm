package PDBx::Raw;

use strict;
use warnings;

use Carp;
use List::Util qw( uniq );

# ------------------------- Constructors/Destructors -------------------------- #

sub new
{
    my ( $class, $pdbx_file, $data_identifier ) = @_;
    my $self = {};
    if( defined $pdbx_file && defined $data_identifier ) {
        $self = _obtain_pdbx_data( $pdbx_file, $data_identifier );
        return bless $self, $class;
    }
    return bless $self, $class;
}

# ----------------------------- Setters/Getters ------------------------------- #

sub get_category
{
    my ( $self, $category ) = @_;
    if( exists $self->{$category} ) {
        return $self->{$category};
    } else {
        return {};
    }
}

sub set_category
{
    my ( $self, $category, $data ) = @_;
    $self->{$category} = $data;
}

# --------------------------------- Methods ----------------------------------- #

sub read
{
    my ( $self, $pdbx_file, $data_identifier ) = @_;
    my $pdbx_data = _obtain_pdbx_data( $pdbx_file, $data_identifier );
    for my $category ( keys %{ $pdbx_data } ) {
        $self->{$category} = $pdbx_data->{$category};
    }
}

sub raw
{
    my ( $self ) = @_;

    my $pdbx_indexed = $self;

    my %pdbx_raw = ();
    for my $category ( keys %{ $pdbx_indexed } ) {
        my $current_pdbx_indexed = $pdbx_indexed->{$category}{'data'};
        my @category_attributes =
            sort { $a cmp $b }
            uniq
            map { keys %{ $current_pdbx_indexed->{$_} } }
            keys %{ $current_pdbx_indexed };

        $pdbx_raw{$category}{'metadata'}{'attributes'} = \@category_attributes;
        $pdbx_raw{$category}{'metadata'}{'is_loop'} =
            $pdbx_indexed->{$category}{'metadata'}{'is_loop'};

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
    }

    return \%pdbx_raw;
}

#
# Converts atom site data structure to PDBx.
# Input:
#     $options->{data_name} - data name of the PDBx;
#     $options->{attributes} - attribute list that should be included in the
#     output;
#     $options->{add_attributes} - add list of attributes to existing data
#     structure;
#     $options->{fh} - file handler.
# Output:
#     PDBx STDOUT.
#

sub to_pdbx
{
    my ( $self, $options ) = @_;
    my ( $data_name, $attributes, $add_attributes, $fh ) = (
        $options->{'data_name'},
        $options->{'attributes'},
        $options->{'add_attributes'},
        $options->{'fh'},
    );

    $data_name //= '';
    $attributes //= { '_atom_site' =>
                          [ 'group_PDB',
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
                            'pdbx_PDB_model_num', ] };
    $add_attributes //= { '_atom_site' => [ '_[local]_selection_state',
                                            '_[local]_selection_group' ] };
    $fh //= \*STDOUT;

    print {$fh} "data_${data_name}\n#\n";

    my $pdbx_data = $self;

    if( %{ $pdbx_data } ) {
        for my $category  ( sort { $a cmp $b } keys %{ $pdbx_data } ) {
            if( $pdbx_data->{$category}{'metadata'}{'is_loop'} ) {
                print {$fh} "loop_\n";
                foreach( @{ $pdbx_data->{$category}{'metadata'}{'attributes'} } ) {
                    print {$fh} "$category.$_\n";
                }
                my $attribute_array_length =
                    $#{ $pdbx_data->{$category}{'metadata'}{'attributes'} };
                my $data_array_length =
                    $#{ $pdbx_data->{$category}{'data'} };
                for( my $i = 0;
                     $i <= $data_array_length;
                     $i += $attribute_array_length + 1 ){
                    my @current_data_list = ();
                    for my $data_value ( @{ $pdbx_data->{$category}{'data'} }
                                          [$i..$i+$attribute_array_length] ) {
                        if( defined $data_value ) {
                            push @current_data_list, $data_value;
                        } else {
                            push @current_data_list, '?';
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
    my ( $self, $options ) = @_;
    my ( $category, $attributes ) = (
        $options->{'category'}, $options->{'attributes'}
    );

    $category //= (sort keys %{ $self })[0]; # TODO: assigning through list
                                             # might loose the information.
    my $pdbx_raw = $self->{$category};

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

# --------------------------------- Static  ----------------------------------- #

#
# Obtains pdbx data for a specified categories or items.
# Input:
#     $pdbx_file - PDBx file path;
#     $data_identifier - list of categories or items;
# Output:
#     %pdbx_data - data structure for pdbx data.
#

sub _obtain_pdbx_data
{
    my ( $pdbx_file, $data_identifier ) = @_;
    my %pdbx_data = ();

    if( defined $data_identifier && @{ $data_identifier } ) {
        %pdbx_data = (
            %pdbx_data,
            %{ _obtain_pdbx_line( $pdbx_file, $data_identifier ) } );
        %pdbx_data = (
            %pdbx_data,
            %{ _obtain_pdbx_loop( $pdbx_file, $data_identifier,
                                  { 'ignore_missing_categories' => 1 } ) } );
    }

    return \%pdbx_data;
}

#
# Obtains pdbx lines for a specified items.
# Input:
#     $pdbx_file - PDBx file path;
#     $items - list of desired items.
# Output:
#     %pdbx_line_data - hash of item values.
#

sub _obtain_pdbx_line
{
    my ( $pdbx_file, $items ) = @_;

    my %pdbx_line_data;
    my %current_line_data;
    my $item_regexp = join q{|}, @{ $items };
    $item_regexp =~ s/\[/\\[/g;
    $item_regexp =~ s/\]/\\]/g;

    local $/ = '';
    local @ARGV = ( $pdbx_file );
    while( <> ) {
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
    }

    return \%pdbx_line_data;
}

#
# Obtains pdbx loops for a specified categories.
# Input:
#     $pdbx_file - PDBx file path;
#     $categories - list of specified categories.
#     $options->{'read_until_end'} - reads whole pdbx file or stdin.
# Output:
#     %pdbx_loop_data - data structure for loop data or list of data structure.
#

sub _obtain_pdbx_loop
{
    my ( $pdbx_file, $categories, $options ) = @_;
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

    local @ARGV = ( $pdbx_file );
    my $line_counter = 0;

    while( <> ) {
        if( /^data_/ || ! @categories ) {
            push @categories, [];
            push @attributes, [];
            push @data, [];
        } elsif( /($category_regexp)[.](.+)\n$/x ) {
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
                    if( $pdbx_file eq '-' ) {
                        warn "'$searched_category' data was not found in " .
                             "STDIN.\n";
                    } else {
                        warn "'$searched_category' data was not found in " .
                             "'$pdbx_file'.\n";
                    }
                }
            }
        }
    }

    warn "$pdbx_file - is empty.\n"
        if $line_counter == 0 && $pdbx_file ne '-' && ! $ignore_missing_categories;
    warn "STDIN - is empty.\n"
        if $line_counter == 0 && $pdbx_file eq '-' && ! $ignore_missing_categories;

    # Generates hash from three lists.
    my @pdbx_loop_data;
    for( my $i = 0; $i <= $#categories; $i++ ) {
        my %pdbx_loop_data;
        for( my $j = 0; $j <= $#{ $categories[-1] }; $j++ ) {
            $pdbx_loop_data{$categories[$i][$j]}{'metadata'}{'is_loop'} = 1;
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

1;
