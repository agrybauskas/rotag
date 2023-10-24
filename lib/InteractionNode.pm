package InteractionNode;

use strict;
use warnings;

sub new
{
    my ( $class, $args ) = @_;
    my $self = {
        "residue_id" => undef,
        "neighbour_residue_id" => undef,
    };
    return bless $self, $class;
}

1;
