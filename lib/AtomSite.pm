package AtomSite;

use strict;
use warnings;

use List::MoreUtils qw( any );
use List::Util qw( max );

use PDBxParser qw( pdbx_loop_unique
                   obtain_pdbx_loop
                   to_pdbx );

# ------------------------- Constructor and destructor ------------------------ #

sub new
{
    my ( $class, $options ) = @_;
    my $self = { 'last_atom_id' => 0 };

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
    my ( $self, $pdbx_file ) = @_;

    $self->{'atoms'} =
        pdbx_loop_unique( obtain_pdbx_loop( $pdbx_file, [ '_atom_site' ] ) );

    $self->{'last_atom_id'} = max( keys %{ $self->{'atoms'} } );

    return;
}

#
# Creates atom data structure item.
# Input:
#     $options - hash of all necessary attributes with corresponding values;
# Output:
#     atom data structure item.
#

sub create
{
    my ( $self, $options ) = @_;
    my $atom_id = $options->{'id'};
    $atom_id //= max( keys %{ $self->{'atoms'} } );
    my $type_symbol = $options->{'type_symbol'};
    my $label_atom_id = $options->{'label_atom_id'};
    my $label_alt_id = $options->{'label_alt_id'};
    $label_alt_id //= q{.};
    my $label_comp_id = $options->{'label_comp_id'};
    my $label_asym_id = $options->{'label_asym_id'};
    my $label_entity_id = $options->{'label_entity_id'};
    $label_entity_id //= q{?};
    my $label_seq_id = $options->{'label_seq_id'};
    my $cartn_x = $options->{'cartn_x'};
    my $cartn_y = $options->{'cartn_y'};
    my $cartn_z = $options->{'cartn_z'};
    my $pdbx_model_num = $options->{'pdbx_PDB_model_num'};

    if( ! exists $self->{'atoms'}{$atom_id} ) {
        $self->{'atoms'}{$atom_id}{'group_PDB'} = 'ATOM';
        $self->{'atoms'}{$atom_id}{'id'} = $atom_id;
        $self->{'atoms'}{$atom_id}{'type_symbol'} = $type_symbol;
        $self->{'atoms'}{$atom_id}{'label_atom_id'} = $label_atom_id;
        $self->{'atoms'}{$atom_id}{'label_alt_id'} = $label_alt_id;
        $self->{'atoms'}{$atom_id}{'label_comp_id'} = $label_comp_id;
        $self->{'atoms'}{$atom_id}{'label_asym_id'} = $label_asym_id;
        $self->{'atoms'}{$atom_id}{'label_entity_id'} = $label_entity_id;
        $self->{'atoms'}{$atom_id}{'label_seq_id'} = $label_seq_id;
        $self->{'atoms'}{$atom_id}{'Cartn_x'} = $cartn_x;
        $self->{'atoms'}{$atom_id}{'Cartn_y'} = $cartn_y;
        $self->{'atoms'}{$atom_id}{'Cartn_z'} = $cartn_z;
        $self->{'atoms'}{$atom_id}{'pdbx_PDB_model_num'} = $pdbx_model_num;
    } else {
        die 'Specified atom id is already present in the atom site';
    }

    return;
}

#
# Appends atom data structure to the current object.
# Input:
#     $atom_sites - list of appendable atom data structures.
#     $action - flags: 'r' - renumber, 'o' - overwrite.
# Output:
#     none - atom sites are appended to the structure.
#

sub append
{
    my ( $self ) = shift;
    my ( $atom_sites, $action ) = @_;
    $action //= '';

    for my $atom_site ( @{ $atom_sites } ) {
        if( $action eq 'r' ) {
            for my $atom_id ( sort keys %{ $atom_site } ) {
                $self->{'atoms'}{$self->{'last_atom_id'}+1} =
                    $atom_site->{$atom_id};
                $self->{'atoms'}{$self->{'last_atom_id'}+1}{'id'} =
                    $self->{'last_atom_id'}+1;
                $self->{'last_atom_id'}++;
            }
        } else {
            for my $atom_id ( sort keys %{ $atom_site } ) {
                if( ! exists $self->{'atoms'}{$atom_id} ) {
                    $self->{'atoms'}{$atom_id} = $atom_site->{$atom_id};
                } else {
                    die "Atom with id $atom_id already exists";
                }
            }
        }
    }

    $self->{'last_atom_id'} = max( keys %{ $self->{'atoms'} } );

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
#     structure.
# Output:
#     \%filtered_atoms- filtered atom data structure;
#

sub filter
{
    my ( $self, $options ) = @_;
    my ( $include, $exclude ) =
        ( $options->{'include'}, $options->{'exclude'} );

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

    # Return object handle.
    my $filtered_atom_site = AtomSite->new();
    $filtered_atom_site->append( [ \%filtered_atoms ] );

    return $filtered_atom_site;
}

#
# Converts atom site data structure to PDBx.
# Input:
#     $options->{'atoms'} - atom site data structure.
# Output:
#     STDERR - PDBx format file.
#

sub pdbx
{
    my ( $self ) = @_;

    return to_pdbx( { 'atom_site' => $self->{'atoms'}, 'data_name' => '' } );
}

1;
