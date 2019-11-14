package PredictSidechains;

use strict;
use warnings;

use Exporter qw( import );
our @EXPORT_OK = qw( sidechain_positions );

use Graph;
use GraphViz;

use PDBxParser qw( filter_new );
use Grid qw( grid_box
             identify_neighbour_cells );

use Version qw( $VERSION );

our $VERSION = $VERSION;

sub sidechain_positions
{
    my ( $args ) = @_;
    my ( $parameters, $atom_site ) =
        ( $args->{'parameters'}, $args->{'atom_site'} );

    my $edge_length_interaction =
        $parameters->{'_[local]_constants'}{'edge_length_interaction'};

    # Chooses only CA atoms, because from them, boundary interaction conditions
    # are measured.
    my $atom_site_cas = filter_new( $atom_site,
                                    { 'include' =>
                                          { 'label_atom_id' => [ 'CA' ] } } );
    my ( $grid_box_cas ) =
        grid_box( $parameters, $atom_site_cas, $edge_length_interaction );

}

# ----------------------------------------------------------------------------- #

1;
