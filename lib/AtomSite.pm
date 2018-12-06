package AtomSite;

use strict;
use warnings;

require Exporter;

our @ISA = qw( Exporter );
our @EXPORT = qw( extract
                  filter
                  mark_selection
                  pdbx );

use Carp;
use List::MoreUtils qw( any
                        uniq );
use List::Util qw( max );

use Atom;
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

# --------------------------------- Methods ----------------------------------- #

#
# From PDBx file, obtains data only from _atom_site category and outputs special
# atom data structure that represents atom data.
# Input:
#     $pdbx_file - PDBx file.
# Output:
#     Atom data structure is created and stored in $self->{'atoms'}.
#

sub open
{
    my ( $self, $pdbx_file ) = @_;

    my $atom_site =
        pdbx_loop_unique( obtain_pdbx_loop( $pdbx_file, [ '_atom_site' ] ) );

    for my $atom_id ( keys %{ $atom_site } ) {
        $self->{'_atoms'}{$atom_id} = Atom->new( $atom_site->{$atom_id} );
    }

    return;
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
# Adds atoms to atom site.
# Input:
#     $atoms - list of atom data structure.
# Output:
#     adds atoms to the atom site.
#

sub add
{
    my ( $self, $atoms ) = @_;

    for my $atom ( @{ $atoms } ) {
        if( ref $atom ne 'Atom' ) {
            confess "not 'Atom' struct is being added";
        }

        my $atom_id = $atom->{'id'};
        if( ! exists $self->{'_atoms'}{$atom_id} ) {
            $self->{'_atoms'}{$atom_id} = $atom;
        } else {
            confess "atom with id $atom_id already exists";
        }
    }

    return;
}

#
# Appends atom data structure to the current object.
# Input:
#     $atom_sites - list of appendable atom site data structures.
# Output:
#     none - atom sites are appended to the structure.
#

sub append
{
    my ( $self, $atom_sites ) = @_;

    for my $atom_site ( @{ $atom_sites } ) {
        if( ref $atom_site ne 'AtomSite' ) {
            confess "not 'AtomSite' object is being appended";
        }

        for my $atom_id ( keys %{ $atom_site->{'_atoms'} } ) {
            if( ! exists $self->{'_atoms'}{$atom_id} ) {
                $self->{'_atoms'}{$atom_id}=$atom_site->{'_atoms'}{$atom_id};
            } else {
                confess "atom with id $atom_id already exists";
            }
        }
    }

    return;
}

#
# Adds T (target), S (surrounding or selected), I (ignored) char to
# '[local]_selection_state' attribute data.
# Input:
#     $atom_site - atom site data structure;
#     $options{'target'} - list of ids of the target atom;
#     $options{'select'} - list of ids of the selected atom;
#     $options{'ignore'} - list of ids of the ignored atom.
# Output:
#     adds markers to specified attribute field.
#

sub mark_selection
{
    my ( $atom_site, $options ) = @_;

    my ( $target_atom_ids, $selected_atom_ids ) =
        ( $options->{'target'}, $options->{'select'}, );

    for my $atom_id ( keys %{ $atom_site->{'_atoms'} } ) {
        $atom_site->{'_atoms'}{$atom_id}{'[local]_selection_state'} = 'I';
    }

    for my $selected_atom_id ( @{ $selected_atom_ids } ) {
        $atom_site->{'_atoms'}{$selected_atom_id}{'[local]_selection_state'}='S';
    }

    for my $target_atom_id ( @{ $target_atom_ids } ) {
        $atom_site->{'_atoms'}{$target_atom_id}{'[local]_selection_state'} = 'T';
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
    my ( $atom_site, $options ) = @_;
    my ( $include, $exclude, $return_data, $group_id, $selection_state ) =
        ( $options->{'include'}, $options->{'exclude'},
          $options->{'return_data'}, $options->{'group_id'},
          $options->{'selection_state'} );

    if( ! defined $atom_site->{'_atoms'} ) {
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

    for my $atom_id ( keys %{ $atom_site->{'_atoms'} } ) {
        my $keep_atom = 1;
        for my $attribute ( keys %{ $atom_site->{'_atoms'}{$atom_id} } ) {
            my $value = $atom_site->{'_atoms'}{$atom_id}{$attribute};
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
            $atom_site->{'_atoms'}{$atom_id}{'[local]_selection_group'} =
                $group_id;
        }

        if( defined $selection_state ) {
            $atom_site->{'_atoms'}{$atom_id}{'[local]_selection_state'} =
                $selection_state;
        }

        $filtered_atoms{$atom_id} = $atom_site->{'_atoms'}{$atom_id};
    }

    # Return object handle or atom ids depending on the flag.
    if( defined $return_data && $return_data eq 'id' ) {
        return [ keys %filtered_atoms ];
    } elsif( defined $return_data ) {
        return [ map { $atom_site->{'_atoms'}{$_}{$return_data} }
                 keys %filtered_atoms ];
    } else {
        return bless { '_atoms' => \%filtered_atoms }, 'AtomSite';
    }
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

        # Simply iterates through $atom_site->{'_atoms'} keys and extracts data
        # using data specifier and is asigned to atom id.
        for my $atom_id ( sort { $a <=> $b } keys %{ $atom_site->{'_atoms'} } ) {
            $atom_data_with_id{$atom_id} =
                [ map { $atom_site->{'_atoms'}->{$atom_id}{$_} } @{ $data } ];
        }

        return \%atom_data_with_id;
    } else {
        for my $atom_id ( sort { $a <=> $b } keys %{ $atom_site->{'_atoms'} } ) {
            if( defined $is_list && $is_list ) {
                push @atom_data,
                    map { $atom_site->{'_atoms'}->{$atom_id}{$_} } @{ $data };
            } else {
                push @atom_data,
                    [ map { $atom_site->{'_atoms'}->{$atom_id}{$_} } @{ $data } ];
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
    my ( $atom_site, $options ) = @_;
    my ( $data_name ) = $options->{'data_name'};
    $data_name = '';

    return to_pdbx( { 'atom_site' => $atom_site->{'_atoms'},
                      'data_name' => $data_name } );
}

1;
