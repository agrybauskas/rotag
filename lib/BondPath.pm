package BondPath;

use strict;
use warnings;

use AtomProperties qw( sort_atom_ids_by_name
                       sort_by_unique_residue_key );
use PDBxParser qw( filter_new );

# ------------------------- Constructors/Destructors -------------------------- #

sub new
{
    my ( $class, $args ) = @_;
    my ( $atom_site, $start_atom_ids, $start_atom_names, $include_hetatoms,
         $ignore_connections, $include_visited ) =
        ( $args->{'atom_site'},
          $args->{'start_atom_ids'},
          $args->{'start_atom_names'},
          $args->{'include_hetatoms'},
          $args->{'ignore_connections'},
          $args->{'include_visited'} );

    my $self = {};

    # By default, N is starting atom for main-chain calculations. XA is added
    # for debugging and test purposes.
    $start_atom_names //= [ 'N', 'XA' ];
    $start_atom_ids //=
        filter_new( $atom_site,
                    { 'include' => { 'label_atom_id' => $start_atom_names },
                      'return_data' => 'id' } );
    $include_hetatoms //= 0;
    $ignore_connections //= {};
    $include_visited //= 0;

    my %visited_atom_ids = (); # Contains visited atom order.
    my @next_atom_ids =
        grep { $atom_site->{$_}{'group_PDB'} ne 'HETATM' }
        grep { defined $_ }
             ( shift @{ sort_by_unique_residue_key( $start_atom_ids,
                                                    $atom_site ) } );

    my %bond_paths = ();
    my $atom_order_idx = 1;

    # Exists if there are no atoms that is not already visited.
    while( @next_atom_ids ) {
        my ( $atom_id ) = pop @next_atom_ids;
        my $atom_name = $atom_site->{$atom_id}{'label_atom_id'};

        next if $visited_atom_ids{$atom_id};
        $visited_atom_ids{$atom_id} = $atom_order_idx;
        $self->{'atom_order'}{$atom_id} = $atom_order_idx;

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
            my $sorted_neighbour_atom_name =
                $atom_site->{$sorted_neighbour_atom_id}{'label_atom_id'};

            next if $visited_atom_ids{$sorted_neighbour_atom_id};

            next if exists $ignore_connections->{'label_atom_id'}{$atom_name} &&
                exists $ignore_connections->{'label_atom_id'}
                                            {$atom_name}
                                            {$sorted_neighbour_atom_name} &&
                $ignore_connections->{'label_atom_id'}
                                     {$atom_name}
                                     {$sorted_neighbour_atom_name};

            unshift @next_atom_ids, $sorted_neighbour_atom_id;
        }

        $atom_order_idx++;
    }

    # Prepares connection information.
    # HACK: needs to be optimized and moved with previous code if possible.
    for my $atom_id ( sort { $self->{'atom_order'}{$a} <=>
                             $self->{'atom_order'}{$b} }
                      keys %{ $self->{'atom_order'} } ) {

        my $atom_name = $atom_site->{$atom_id}{'label_atom_id'};

        my @atom_connections =
            map { [ $_, 'connections' ] }
               @{ $atom_site->{$atom_id}{'connections'} };
        if( $include_hetatoms &&
            defined $atom_site->{$atom_id}{'connections_hetatom'} ){
            push @atom_connections,
                map { [ $_, 'connections_hetatom' ] }
                   @{ $atom_site->{$atom_id}{'connections_hetatom'} };
        }

        next if ! @atom_connections;

        for my $atom_connection ( @atom_connections ) {
            my $neighbour_atom_id = $atom_connection->[0];
            my $connection_type = $atom_connection->[1];

            next if ! exists $atom_site->{$neighbour_atom_id};

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

            if( ! exists $self->{'to'}{$neighbour_atom_id} ||
                $self->{'atom_order'}{$atom_id} <
                $self->{'atom_order'}{$self->{'to'}{$neighbour_atom_id}} ) {
                $self->{'to'}{$neighbour_atom_id} = $atom_id;
            }

            $self->{'connection_type'}{$atom_id}{$neighbour_atom_id} =
                $connection_type;
            $self->{'connection_type'}{$neighbour_atom_id}{$atom_id} =
                $connection_type;
        }
    }

    return bless $self, $class;
}

# ----------------------------- Setters/Getters ------------------------------- #

sub get_atom_order
{
    my ( $self ) = @_;
    return [ sort { $self->{'atom_order'}{$a} <=> $self->{'atom_order'}{$b} }
             keys %{ $self->{'atom_order'} } ];
}

sub get_atom_id_to
{
    my ( $self, $atom_id ) = @_;
    return $self->{'to'}{$atom_id}
}

sub get_connection_type
{
    my ( $self, $atom_id_from, $atom_id_to ) = @_;
    return $self->{'connection_type'}{$atom_id_from}{$atom_id_to};
}

# --------------------------------- Methods ----------------------------------- #

1;
