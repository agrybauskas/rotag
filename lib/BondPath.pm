package BondPath;

use strict;
use warnings;

use AtomProperties qw( sort_atom_ids_by_name );
use PDBxParser qw( filter_new );

# ------------------------- Constructors/Destructors -------------------------- #

sub new
{
    my ( $class, $args ) = @_;
    my ( $atom_site, $start_atom_ids, $include_hetatoms ) =
        ( $args->{'atom_site'},
          $args->{'start_atom_ids'},
          $args->{'include_hetatoms'} );

    my $self = {};

    # By default, N is starting atom for main-chain calculations. XA is added
    # for debugging and test purposes.
    # TODO: think about more general rule -- how to choose the starting atom.
    $start_atom_ids //=
        filter_new( $atom_site,
                    { 'include' => { 'label_atom_id' => [ 'N', 'XA' ] },
                      'return_data' => 'id' } );

    my %visited_atom_ids = ();
    my @next_atom_ids = ( grep { defined $_ } @{ $start_atom_ids } );

    my %bond_paths = ();
    my $atom_order_idx = 1;

    # Exists if there are no atoms that is not already visited.
    while( @next_atom_ids ) {
        my ( $atom_id ) = pop @next_atom_ids;

        next if $visited_atom_ids{$atom_id};
        $visited_atom_ids{$atom_id} = 1;

        # Marks neighbouring atoms.
        my @neighbour_atom_ids = ();
        if( defined $atom_site->{$atom_id}{'connections'} ) {
            push @neighbour_atom_ids,
                grep { defined $atom_site->{$_} }
                    @{ $atom_site->{$atom_id}{'connections'} };
        }
        if( $include_hetatoms &&
            defined $atom_site->{$atom_id}{'connections_hetatom'} ) {
            push @neighbour_atom_ids,
                grep { defined $atom_site->{$_} }
                    @{ $atom_site->{$atom_id}{'connections_hetatom'} };
        }

        my @sorted_neighbour_atom_ids =
            @{ sort_atom_ids_by_name( \@neighbour_atom_ids, $atom_site ) };

        for( my $i = 0; $i <= $#sorted_neighbour_atom_ids; $i++ ) {
            my $sorted_neighbour_atom_id = $sorted_neighbour_atom_ids[$i];

            next if $visited_atom_ids{$sorted_neighbour_atom_id};

            if( ! exists $bond_paths{$sorted_neighbour_atom_id} ) {
                $self->{'connections'}{$atom_id}{$sorted_neighbour_atom_id} =
                    { 'order' => $atom_order_idx };
                $self->{'order'}{$atom_order_idx} =
                    [ $atom_id, $sorted_neighbour_atom_id ];
            }

            # Depth-first search.
            if( $i == 0 ) {
                push @next_atom_ids, $sorted_neighbour_atom_id;
            } else {
                unshift @next_atom_ids, $sorted_neighbour_atom_id;
            }

            $atom_order_idx++;
        }
    }

    return bless $self, $class;
}

1;
