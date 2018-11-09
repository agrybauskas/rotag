package AtomSite;

use strict;
use warnings;

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
# data structure that represents atom data.
# Input:
#     $pdbx_file - PDBx file.
# Output:
#     %atom_site - special data structure.
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

sub get_atom_site
{
    my ( $self ) = @_;

    return $self->{'atom_site'};
}

1;
