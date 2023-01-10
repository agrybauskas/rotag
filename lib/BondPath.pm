package BondPath;

use strict;
use warnings;

use AtomProperties qw( sort_atom_ids_by_name );
use PDBxParser qw( filter_new );

# ------------------------- Constructors/Destructors -------------------------- #

sub new
{
    my ( $class, $args ) = @_;
    my ( $atom_site, $start_atom_ids, $include_hetatoms, $ignore_connections,
         $include_visited ) =
        ( $args->{'atom_site'},
          $args->{'start_atom_ids'},
          $args->{'include_hetatoms'},
          $args->{'ignore_connections'},
          $args->{'include_visited'} );

    my $self = {};

    # By default, N is starting atom for main-chain calculations. XA is added
    # for debugging and test purposes.
    # TODO: think about more general rule -- how to choose the starting atom.
    $start_atom_ids //=
        filter_new( $atom_site,
                    { 'include' => { 'label_atom_id' => [ 'N', 'XA' ] },
                      'return_data' => 'id' } );
    $include_hetatoms //= 0;
    $ignore_connections //= {};
    $include_visited //= 0;

    my %visited_atom_ids = (); # Contains visited atom order.
    my @next_atom_ids = ( grep { defined $_ } @{ $start_atom_ids } );

    my %bond_paths = ();
    my $atom_order_idx = 1;

    # Exists if there are no atoms that is not already visited.
    while( @next_atom_ids ) {
        my ( $atom_id ) = pop @next_atom_ids;

        next if $visited_atom_ids{$atom_id};
        $visited_atom_ids{$atom_id} = $atom_order_idx;
        $self->{'atom_order'}{$atom_order_idx} = $atom_id;

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

            # Depth-first search.
            if( $i == 0 ) {
                push @next_atom_ids, $sorted_neighbour_atom_id;
            } else {
                unshift @next_atom_ids, $sorted_neighbour_atom_id;
            }
        }

        $atom_order_idx++;
    }

    # Prepares connection information.
    for my $order ( sort { $a <=> $b } keys %{ $self->{'atom_order'} } ) {
        my $atom_id = $self->{'atom_order'}{$order};
        my $atom_name = $atom_site->{$atom_id}{'label_atom_id'};
        my $atom_connections = $atom_site->{$atom_id}{'connections'};

        next if ! defined $atom_connections;

        for my $neighbour_atom_id ( @{ $atom_connections } ) {
            my $neighbour_atom_name =
                $atom_site->{$neighbour_atom_id}{'label_atom_id'};

            next if exists $ignore_connections->{'label_atom_id'}{$atom_name} &&
                exists $ignore_connections->{'label_atom_id'}
                                            {$atom_name}
                                            {$neighbour_atom_name} &&
                $ignore_connections->{'label_atom_id'}
                                     {$atom_name}
                                     {$neighbour_atom_name};

            next if ! $include_visited &&
                $visited_atom_ids{$atom_id} > $visited_atom_ids{$neighbour_atom_id};

            $self->{'connections'}{$atom_id}{$neighbour_atom_id} =
                $visited_atom_ids{$neighbour_atom_id};
        }
    }

    return bless $self, $class;
}

1;
