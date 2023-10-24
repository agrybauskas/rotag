package InteractionNode;

use strict;
use warnings;

sub new
{
    my ( $class, $args ) = @_;
    my $self = {
        "residue_id" => undef,
        "neighbour_ids" => undef,
        "visited" => 0,
    };
    return bless $self, $class;
}

sub residue_id
{
    my ( $self, $residue_id ) = @_;
    if( ! defined $residue_id ) {
        return $self->{'residue_id'};
    }
    $self->{'residue_id'} = $residue_id;
}

sub set_visited
{
    my ( $self ) = @_;
    $self->{'visited'} = 1;
}

sub is_visited
{
    my ( $self ) = @_;
    return $self->{'visited'};
}

1;
