package AtomSite;

use strict;
use warnings;

use PDBxParser qw( obtain_atom_site );

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

sub open
{
    my ( $self, $pdbx_file ) = @_;

    $self->{'atom_site'} = obtain_atom_site( $pdbx_file );

    return;
}

sub get_atom_site
{
    my ( $self ) = @_;

    return $self->{'atom_site'};
}

1;
