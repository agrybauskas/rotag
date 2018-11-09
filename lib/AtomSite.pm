package AtomSite;

use strict;
use warnings;

use List::MoreUtils qw( any );

use PDBxParser qw( pdbx_loop_unique
                   obtain_pdbx_loop );

# ------------------------- Constructor and destructor ------------------------ #

sub new
{
    my ( $class, $options ) = @_;
    my $self = {};

    return bless $self, $class;
}

sub destroy
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
#     none - atom data structure is created and stored in $self->{'atom_site'}.
#
#     Ex.: { 1 => { 'group_id' => 'ATOM',
#                   'id'       => 1,
#                   ... } }
#

sub open
{
    my ( $self, $pdbx_file ) = @_;

    $self->{'atom_site'} =
        pdbx_loop_unique( obtain_pdbx_loop( $pdbx_file, [ '_atom_site' ] ) );

    return;
}

#
# Returns atom data structure.
# Output:
#     atom data structure.
#

sub get_atom_site
{
    my ( $self ) = @_;

    return $self->{'atom_site'};
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

    my $atom_site = $self->{'atom_site'};

    if( ! defined $self->{'atom_site'} ) {
        die 'No atom were loaded to the AtomSite data structure.';
    }

    # Iterates through each atom in $self->{'atom_site'} and checks if atom
    # specifiers match up.
    my %filtered_atom_site;

    # First, filters atoms that are described in $self->{'include'} specifier.
    if( defined $include ) {
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
                $filtered_atom_site{$atom_id} = $atom_site->{$atom_id};
            }
        }
    } else {
        %filtered_atom_site = %{ $atom_site };
    }

    # Then filters out atoms that are in $self->{'exclude'} specifier.
    if( defined $exclude ) {
        for my $atom_id ( keys %filtered_atom_site ) {
            for my $attribute ( keys %{ $exclude } ) {
                if( exists $atom_site->{$atom_id}{$attribute} &&
                    any { $atom_site->{$atom_id}{$attribute} eq $_ }
                       @{ $exclude->{$attribute} } ) {
                    delete $filtered_atom_site{$atom_id};
                    last;
                }
            }
        }
    }

    return \%filtered_atom_site;
}

1;
