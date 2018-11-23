package AtomSite;

use strict;
use warnings;

use List::MoreUtils qw( any
                        uniq );
use List::Util qw( max );

use BondProperties qw( hybridization );
use ConnectAtoms qw( connect_atoms );
use PDBxParser qw( pdbx_loop_unique
                   obtain_pdbx_loop
                   to_pdbx );

# ------------------------- Constructor and destructor ------------------------ #

sub new
{
    my ( $class, $options ) = @_;

    my $self = { '_atoms' => {} };

    return bless $self, $class;
}

sub DESTROY
{
    my ( $self ) = @_;
}

# --------------------------------- Methods ----------------------------------- #

#
# From PDBx file, obtains data only from _atom_site category and outputs special
# atom data structure that represents atom data.
# Input:
#     $pdbx_file - PDBx file.
# Output:
#     none - atom data structure is created and stored in $self->{'atoms'}.
#
#     Ex.: { 1 => { 'group_id' => 'ATOM',
#                   'id'       => 1,
#                   ... } }
#

sub open
{
    my ( $self, $pdbx_file, $options ) = @_;

    $self->{'_atoms'} =
        pdbx_loop_unique( obtain_pdbx_loop( $pdbx_file, [ '_atom_site' ] ) );

    $self->{'_last_atom_id'} = max( keys %{ $self->{'_atoms'} } );

    return;
}

#
# Creates atom data structure item.
# Input:
#     $atom_data - hash of all necessary attributes with corresponding values;
# Output:
#     %atom_site - atom site data structure.
#

sub create
{
    my ( $self, $atom_data ) = @_;
    my $atom_id = $atom_data->{'id'};
    my $type_symbol = $atom_data->{'type_symbol'};
    my $label_atom_id = $atom_data->{'label_atom_id'};
    my $label_alt_id = $atom_data->{'label_alt_id'};
    $label_alt_id //= q{.};
    my $label_comp_id = $atom_data->{'label_comp_id'};
    my $label_asym_id = $atom_data->{'label_asym_id'};
    my $label_entity_id = $atom_data->{'label_entity_id'};
    $label_entity_id //= q{?};
    my $label_seq_id = $atom_data->{'label_seq_id'};
    my $cartn_x = $atom_data->{'Cartn_x'};
    my $cartn_y = $atom_data->{'Cartn_y'};
    my $cartn_z = $atom_data->{'Cartn_z'};
    my $pdbx_model_num = $atom_data->{'pdbx_PDB_model_num'};

    my %atom_site = ();
    if( ! exists $atom_site{$atom_id} ) {
        $atom_site{$atom_id}{'group_PDB'} = 'ATOM';
        $atom_site{$atom_id}{'id'} = $atom_id;
        $atom_site{$atom_id}{'type_symbol'} = $type_symbol;
        $atom_site{$atom_id}{'label_atom_id'} = $label_atom_id;
        $atom_site{$atom_id}{'label_alt_id'} = $label_alt_id;
        $atom_site{$atom_id}{'label_comp_id'} = $label_comp_id;
        $atom_site{$atom_id}{'label_asym_id'} = $label_asym_id;
        $atom_site{$atom_id}{'label_entity_id'} = $label_entity_id;
        $atom_site{$atom_id}{'label_seq_id'} = $label_seq_id;
        $atom_site{$atom_id}{'Cartn_x'} = $cartn_x;
        $atom_site{$atom_id}{'Cartn_y'} = $cartn_y;
        $atom_site{$atom_id}{'Cartn_z'} = $cartn_z;
        $atom_site{$atom_id}{'pdbx_PDB_model_num'} = $pdbx_model_num;
    } else {
        die 'Specified atom id is already present in the atom site';
    }

    return \%atom_site;
}

#
# Updates information about bonds in AtomSite object.
# Output:
#     updates/inserts 'connections', 'hybridization' attributes in atom site.
#

sub update_bonds
{
    my ( $self ) = @_;

    connect_atoms( $self->{'_atoms'} );
    hybridization( $self->{'_atoms'} );

    return;
}

#
# Appends atom data structure to the current object.
# Input:
#     $atom_sites - list of appendable atom data structures.
#     $renumber - renumbers the ids.
# Output:
#     none - atom sites are appended to the structure.
#

sub append
{
    my ( $self, $atom_sites ) = @_;

    for my $atom_site ( @{ $atom_sites } ) {
        for my $atom_id ( sort keys %{ $atom_site->{'_atoms'} } ) {
            if( ! exists $self->{'_atoms'}{$atom_id} ) {
                $self->{'_atoms'}{$atom_id} = $atom_site->{'_atoms'}{$atom_id};
            } else {
                die "Atom with id $atom_id already exists";
            }
        }
    }

    return;
}

#
# Creates an index table that stores atom ids according to the attribute data.
# Input:
#     $atom_site - atom site data structure.
# Output:
#     %index_table - creates an index table.
# Ex.:
#     { 'label_atom_id' => { 'CA' => [ 1, 4, 7 ],
#                            'CB' => [ 2, 5, 8 ] } }
#

sub index
{
    my ( $self, $atom_site ) = @_;

    my %index_table;

    for my $atom_id ( keys %{ $atom_site->{'_atoms'} } ) {
        my $atom = $atom_site->{'_atoms'}{$atom_id};

        for my $attribute ( keys %{ $atom } ) {
            my $value = $atom->{$attribute};
            next if ref $value;
            if( exists $index_table{"$attribute"} &&
                exists $index_table{"$attribute"}{"$value"} ) {
                push @{ $index_table{"$attribute"}{"$value"} }, $atom_id;
            } else {
                $index_table{"$attribute"}{"$value"} = [ $atom_id ];
            }
        }
    }

    return \%index_table;
}

#
# Adds T (target), S (surrounding or selected), I (ignored) char to
# '[local]_selection_state' attribute data.
# Input:
#     $options{'target'} - list of ids of the target atom;
#     $options{'select'} - list of ids of the selected atom;
#     $options{'ignore'} - list of ids of the ignored atom.
# Output:
#     adds markers to specified attribute field.
#

sub mark_selection
{
    my ( $self, $options ) = @_;

    my ( $target_atom_ids, $selected_atom_ids ) =
        ( $options->{'target'}, $options->{'select'}, );

    for my $atom_id ( keys %{ $self->{'atoms'} } ) {
        if( any { $atom_id eq $_  } @{ $target_atom_ids } ) {
            $self->{'atoms'}{$atom_id}{'[local]_selection_state'} = 'T';
        } elsif( any { $atom_id eq $_  } @{ $selected_atom_ids } ) {
            $self->{'atoms'}{$atom_id}{'[local]_selection_state'} = 'S';
        } else {
            $self->{'atoms'}{$atom_id}{'[local]_selection_state'} = 'I';
        }
    }

    return;
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

sub filter
{
    my ( $self, $options ) = @_;
    my ( $include, $exclude, $return_ids, $group_id ) =
        ( $options->{'include'},    $options->{'exclude'},
          $options->{'return_ids'}, $options->{'group_id'} );

    $return_ids //= 0;

    my $atoms = $self->{'atoms'};

    if( ! defined $self->{'atoms'} ) {
        die 'No atom were loaded to the AtomSite data structure';
    }

    # Iterates through each atom in $self->{'atoms'} and checks if atom
    # specifiers match up.
    my %filtered_atoms;

    # First, filters atoms that are described in $self->{'include'} specifier.
    if( defined $include ) {
        for my $atom_id ( keys %{ $atoms } ) {
            my $match_counter = 0; # Tracks if all matches occured.
            for my $attribute ( keys %{ $include } ) {
                if( exists $atoms->{$atom_id}{$attribute} &&
                    any { $atoms->{$atom_id}{$attribute} eq $_ }
                       @{ $include->{$attribute} } ) {
                    $match_counter += 1;
                } else {
                    last; # Terminates early if no match is found in specifier.
                }
            }

            if( $match_counter == scalar keys %{ $include } ) {
                $filtered_atoms{$atom_id} = $atoms->{$atom_id};

                # Assigns group id.
                if( defined $group_id ) {
                    $filtered_atoms{$atom_id}{'[local]_selection_group'} =
                        $group_id;
                }
            }
        }
    } else {
        %filtered_atoms = %{ $atoms };
    }

    # Then filters out atoms that are in $self->{'exclude'} specifier.
    if( defined $exclude ) {
        for my $atom_id ( keys %filtered_atoms ) {
            for my $attribute ( keys %{ $exclude } ) {
                if( exists $atoms->{$atom_id}{$attribute} &&
                    any { $atoms->{$atom_id}{$attribute} eq $_ }
                       @{ $exclude->{$attribute} } ) {
                    delete $filtered_atoms{$atom_id};
                    last;
                }
            }
        }
    }

    # Return object handle or atom ids depending on the flag.
    if( $return_ids ) {
        return [ keys %filtered_atoms ]
    } else {
        my $filtered_atom_site = AtomSite->new();
        $filtered_atom_site->append( [ \%filtered_atoms ] );

        return $filtered_atom_site;
    }
}

#
# Extracts information from atom data structure.
# Input:
#     $options{'data'} -
#     $options{'data_with_id'} -
#     $options{'is_list'} -
# Output:
#

sub extract
{
    my ( $self, $options ) = @_;
    my ( $data, $data_with_id, $is_list ) =
        ( $options->{'data'}, $options->{'data_with_id'}, $options->{'is_list'});

    my @atom_data;
    my $atom_site = $self->{'atoms'};

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
# Converts atom site data structure to PDBx.
# Input:
#     $self->{'_atoms'} - atom site data structure.
#     $options->{'data_name'}
# Output:
#     STDERR - PDBx format file.
#

sub pdbx
{
    my ( $self, $options ) = @_;
    my ( $data_name ) = $options->{'data_name'};
    $data_name = '';

    return to_pdbx( { 'atom_site' => $self->{'_atoms'},
                      'data_name' => $data_name } );
}

1;
