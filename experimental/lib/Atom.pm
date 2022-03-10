package Atom;

use strict;
use warnings;

use Carp;

sub new
{
    my ( $class, $atom_data ) = @_;

    # Validates atom data structure.
    my @mandatory_attributes =
        qw( id type_symbol label_atom_id Cartn_x Cartn_y Cartn_z );
    if( ref $atom_data eq 'HASH' ) {
        for my $mandatory_attribute ( @mandatory_attributes ) {
            if( ! exists $atom_data->{$mandatory_attribute} ) {
                confess "'$mandatory_attribute' key is absent in the data " .
                        "structure";
            }
        }
    } else {
        confess "function argument does not have the valid data structure"
    }

    my $self = $atom_data;

    return bless $self, $class;
}

1;
