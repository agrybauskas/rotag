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

    my ( $grid_box_cas, $grid_ca_atom_pos ) =
        grid_box( $parameters, $atom_site_cas, $edge_length_interaction,
                  extract( $atom_site_cas,
                           { 'data' => [ 'id' ], 'is_list' => 1 } ) );

    my $neighbouring_cells_cas =
        identify_neighbour_cells( $grid_box_cas, $grid_ca_atom_pos );

    my %residue_pairs = ();
    for my $grid_id ( keys %{ $grid_ca_atom_pos } ) {
        for my $atom_id ( @{ $grid_ca_atom_pos->{$grid_id} } ) {
            my $unique_residue_key =
                unique_residue_key( $atom_site_cas->{$atom_id} );
            for my $neighbour_atom_id ( @{ $neighbouring_cells_cas->{$grid_id} } ) {
                my $neighbour_unique_residue_key =
                    unique_residue_key( $atom_site_cas->{$neighbour_atom_id} );

                next if $unique_residue_key eq $neighbour_unique_residue_key ||
                    ( defined $residue_pairs{$unique_residue_key} &&
                      defined $residue_pairs{$unique_residue_key}
                                            {$neighbour_unique_residue_key} &&
                      $residue_pairs{$unique_residue_key}
                                    {$neighbour_unique_residue_key} );

                $residue_pairs{$unique_residue_key}
                              {$neighbour_unique_residue_key} = 1;
            }
        }
    }

    # Residues sorted by rotamer number.
    my @sorted_unique_residue_keys =
        sort { scalar( keys %{ $rotamer_look_up_tbls{'unique_residue_key'}
                                                    {$a}{'rotamer_id'} } ) cmp
               scalar( keys %{ $rotamer_look_up_tbls{'unique_residue_key'}
                                                    {$b}{'rotamer_id'} } ) }
        keys %{ $rotamer_look_up_tbls{'unique_residue_key'} };

    # Least-rotamer search in side-chain pairs.
    my $combination_id = 0;
    my @all_residues = @sorted_unique_residue_keys;
    my @next_residues = ( $all_residues[0] );
    while( @next_residues ) {
        @next_residues = (); # Reseting.
    }
}

1;
