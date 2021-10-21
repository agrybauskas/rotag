package PredictSidechains;

use strict;
use warnings;

use Exporter qw( import );
our @EXPORT_OK = qw( predict_sidechains );

use Clone qw( clone );

use Measure qw( energy );
use PDBxParser qw( extract
                   filter_new
                   unique_residue_key
                   to_pdbx );
use PseudoAtoms qw( replace_with_rotamer );
use Grid qw( grid_box
             identify_neighbour_cells );

use Version qw( $VERSION );

our $VERSION = $VERSION;

# ----------------------------- Simple functions ------------------------------ #

sub predict_sidechains
{
    my ( $args ) = @_;

    my ( $atom_site, $rotamer_energies, $rotamer_angles, $parameters ) = (
        $args->{'atom_site'}, $args->{'rotamer_energies'},
        $args->{'rotamer_angles'}, $args->{'parameters'}
    );

    # Error messages for missing arguments.
    if( ! defined $atom_site ) {
        die "atom site is not supplied.\n";
    }
    if( ! defined $rotamer_angles ) {
        die "rotamer angles are not supplied by '_[local]_rotamer_angle' tag.\n";
    }
    if( ! defined $rotamer_energies ) {
        die "rotamer energies are not supplied by '_[local]_rotamer_energy'" .
            " tag.\n";
    }
    if( ! defined $parameters ) {
        die "parameters are not set.\n";
    }

    # Generates look up tables for easier data accesibility.
    my %rotamer_look_up_tbls = ();
    for my $rotamer_angle_id ( keys %{ $rotamer_angles } ) {
        my $rotamer_angle = $rotamer_angles->{$rotamer_angle_id};
        my $rotamer_id = $rotamer_angle->{'rotamer_id'};
        my $unique_residue_key =
            sprintf '%s,%s,%s,%s',
            $rotamer_angle->{'label_seq_id'},
            $rotamer_angle->{'label_asym_id'},
            $rotamer_angle->{'pdbx_PDB_model_num'},
            $rotamer_angle->{'label_alt_id'};

        $rotamer_look_up_tbls{'unique_residue_key'}{$unique_residue_key}
                             {'rotamer_id'}{$rotamer_id} = 1;

        push @{ $rotamer_look_up_tbls{'rotamer_id'}{$rotamer_id}
                                     {'angle_id'} },
            $rotamer_angle_id;
    }

    # Determining interaction grid.
    my $edge_length_interaction =
        $parameters->{'_[local]_constants'}{'edge_length_interaction'};

    # Chooses only CA atoms, because from them, boundary interaction conditions
    # are measured.
    my $atom_site_cas = filter_new( $atom_site,
                                    { 'include' =>
                                          { 'label_atom_id' => [ 'CA' ] } } );

    # TODO: grid box should be built once.
    my ( $grid_box, $grid_atom_pos ) =
        grid_box( $parameters, $atom_site, $edge_length_interaction,
                  extract( $atom_site_cas,
                           { 'data' => [ 'id' ], 'is_list' => 1 } ) );
    my ( $grid_box_cas, $grid_ca_atom_pos ) =
        grid_box( $parameters, $atom_site_cas, $edge_length_interaction,
                  extract( $atom_site_cas,
                           { 'data' => [ 'id' ], 'is_list' => 1 } ) );
    my $neighbouring_cells =
        identify_neighbour_cells( $grid_box, $grid_ca_atom_pos );
    my $neighbouring_cells_cas =
        identify_neighbour_cells( $grid_box_cas, $grid_ca_atom_pos );

    my %grid_ca_atom_pos_by_unique_key = ();
    for my $grid_id ( keys %{ $grid_ca_atom_pos } ) {
        for my $atom_id ( @{ $grid_ca_atom_pos->{$grid_id} } ) {
            my $unique_residue_key =
                unique_residue_key( $atom_site_cas->{$atom_id} );
            if ( ! exists $grid_ca_atom_pos_by_unique_key{$unique_residue_key} ||
                 ! defined $grid_ca_atom_pos_by_unique_key{$unique_residue_key}){
                $grid_ca_atom_pos_by_unique_key{$unique_residue_key} =
                    $grid_id;
                # Adds quick acces in the look up table.
                $rotamer_look_up_tbls{'unique_residue_key'}{$unique_residue_key}
                                     {'grid_id'} = $grid_id;
            }
        }
    }

    # Sort by rotamer number.

    # Least-rotamer search in grid box.
    my $combination_id = 0;
    # TODO: choose better starting point.
    my @all_grid_cells = sort keys %{ $grid_box_cas };
    my @next_grid_cells = ( $all_grid_cells[0] );
    my %visited_grid_cell = ( map { ( $_ => 1 ) } @next_grid_cells );
    while( @next_grid_cells ) {
        my @grid_cells = @next_grid_cells;
        @next_grid_cells = (); # Reseting.
        for my $grid_cell ( @next_grid_cells ) {
            $visited_grid_cell{$grid_cell} = 1;
        }

        # Selecting next grid cells.
    }

    print 'Number of all grid cells: ', scalar @all_grid_cells, "\n";
    print 'Number of Visited grid cells: ',
        scalar( grep { $visited_grid_cell{$_}  } keys %visited_grid_cell ), "\n";
}

1;
