package SidechainInteractions;

use strict;
use warnings;

use Exporter qw( import );

use Clone qw( clone );

use BondProperties qw( hybridization );
use ConnectAtoms qw( connect_atoms );
use Logging qw( info );
use ForceField::Parameters;
use ForceField::Bonded qw( general );
use ForceField::NonBonded qw( general );
use Measure qw( bond_length );
use PDBxParser qw( extract
                   filter_new
                   filter_by_unique_residue_key
                   unique_residue_key
                   to_pdbx );
use PseudoAtoms qw( pairwise_rotamer_energy
                    replace_with_rotamer );
use Grid qw( grid_box
             identify_neighbour_cells );
use SidechainModels qw( rotation_translation );

use Version qw( $VERSION );

our $VERSION = $VERSION;

sub new
{
    my ( $class, $parameters, $atom_site, $rotamer_angles,
         $rotamer_energies ) = @_;

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

    my $self = { 'parameters' => $parameters,
                 'residue_pairs' => undef,
                 'residue_atom_site' => undef,
                 'rotamer_pairs' => undef,
                 'rotamer_atom_site' => undef,
                 'rotamer_angles' => undef,
                 'rotamer_energies' => undef };

    my %residue_to_rotamers = ();
    for my $rotamer_angle_id ( keys %{ $rotamer_angles } ) {
        my $rotamer_angle = $rotamer_angles->{$rotamer_angle_id};
        my $rotamer_id = $rotamer_angle->{'rotamer_id'};
        my $unique_residue_key =
            sprintf '%s,%s,%s,%s',
            $rotamer_angle->{'label_seq_id'},
            $rotamer_angle->{'label_asym_id'},
            $rotamer_angle->{'pdbx_PDB_model_num'},
            $rotamer_angle->{'label_alt_id'};

        $residue_to_rotamers{$unique_residue_key}{$rotamer_id} = 1;
        $self->{'rotamer_angles'}{$rotamer_id}{$rotamer_angle_id} =
            $rotamer_angle;
    }

    # Determining interaction grid.
    my $edge_length_interaction =
        $parameters->{'_[local]_constants'}{'edge_length_interaction'};

    # Chooses only CA atoms, because from them, boundary interaction conditions
    # are measured.
    my $atom_site_cas =
        filter_new( $atom_site,
                    { 'include' =>
                      { 'label_atom_id' => [ 'CA' ] } } );
    my ( $grid_box_cas, $grid_ca_atom_pos ) =
        grid_box( $parameters, $atom_site_cas, $edge_length_interaction,
                  extract( $atom_site_cas,
                           { 'data' => [ 'id' ], 'is_list' => 1 } ) );
    my $neighbouring_cells_cas =
        identify_neighbour_cells( $grid_box_cas, $grid_ca_atom_pos );

    for my $grid_id ( keys %{ $grid_ca_atom_pos } ) {
        for my $atom_id ( @{ $grid_ca_atom_pos->{$grid_id} } ) {
            my $unique_residue_key =
                unique_residue_key( $atom_site->{$atom_id} );

            if( ! exists $self->{'residue_atom_site'}{$unique_residue_key} ) {
                my $residue_site =
                    filter_by_unique_residue_key( $atom_site,
                                                  $unique_residue_key,
                                                  1 );

                connect_atoms( $parameters, $residue_site );
                hybridization( $parameters, $residue_site );
                rotation_translation( $parameters, $residue_site );

                $self->{'residue_atom_site'}{$unique_residue_key} =
                    $residue_site;
            }

            my @rotamer_ids =
                keys %{ $residue_to_rotamers{$unique_residue_key} };

            for my $neighbour_atom_id (@{$neighbouring_cells_cas->{$grid_id}}){
                my $neighbour_unique_residue_key =
                    unique_residue_key( $atom_site->{$neighbour_atom_id} );

                next if $unique_residue_key eq $neighbour_unique_residue_key;

                $self->{'residue_pairs'}{$unique_residue_key}
                                        {$neighbour_unique_residue_key} = 1;
                $self->{'residue_pairs'}{$neighbour_unique_residue_key}
                                        {$unique_residue_key} = 1;

                if( ! exists $self->{'residue_atom_site'}
                                    {$neighbour_unique_residue_key} ) {
                    my $neighbour_residue_site =
                        filter_by_unique_residue_key(
                            $atom_site,
                            $neighbour_unique_residue_key,
                            1
                        );

                    connect_atoms( $parameters, $neighbour_residue_site );
                    hybridization( $parameters, $neighbour_residue_site );
                    rotation_translation( $parameters, $neighbour_residue_site );

                    $self->{'residue_atom_site'}{$neighbour_unique_residue_key}=
                        $neighbour_residue_site;
                }

                my $bond_length = bond_length(
                    [ [ map { $atom_site->{$atom_id}{$_} }
                            ( 'Cartn_x', 'Cartn_y', 'Cartn_z' ) ],
                      [ map { $atom_site->{$neighbour_atom_id}{$_} }
                            ( 'Cartn_x', 'Cartn_y', 'Cartn_z' ) ] ]
                );

                next if $bond_length > $edge_length_interaction;

                my @neighbour_rotamer_ids =
                    keys %{$residue_to_rotamers{$neighbour_unique_residue_key}};

                for my $rotamer_id ( @rotamer_ids ) {
                    for my $neighbour_rotamer_id ( @neighbour_rotamer_ids ) {
                        $self->{'rotamer_pairs'}
                               {$unique_residue_key}
                               {$rotamer_id}
                               {$neighbour_rotamer_id} =
                            1;
                        $self->{'rotamer_pairs'}
                               {$neighbour_unique_residue_key}
                               {$neighbour_rotamer_id}
                               {$rotamer_id} =
                            1;
                    }
                }
            }
        }
    }

    return bless $self, $class;
}

# --------------------------------- Methods ---------------------------------- #

sub predict
{
    my ( $self, $options ) = @_;
    my ( $parameters, $residue_pairs, $rotamer_pairs, $rotamer_angles,
         $residue_atom_site, $rotamer_atom_site ) =
        ( $self->{'parameters'},
          $self->{'residue_pairs'},
          $self->{'rotamer_pairs'},
          $self->{'rotamer_angles'},
          $self->{'residue_atom_site'},
          $self->{'rotamer_atom_site'} );

    my ( $non_bonded_potential, $bonded_potential, $program_called_by,
         $dry_run, $verbose ) =
        ( $options->{'non_bonded_potential'},
          $options->{'bonded_potential'},
          $options->{'program_called_by'},
          $options->{'dry_run'},
          $options->{'verbose'} );

    $dry_run //= 0;
    $verbose //= 0;

    my $cutoff_atom = $parameters->{'_[local]_force_field'}{'cutoff_atom'};

    my %visited_residue_pairs = ();
    my $step_count = 1;
    my @sorted_unique_residue_keys =
        map { $_ }
        sort { scalar( keys %{ $rotamer_pairs->{$a} } ) <=>
               scalar( keys %{ $rotamer_pairs->{$b} } ) }
        keys %{ $residue_pairs };

    while( @sorted_unique_residue_keys ) {
        my $unique_residue_key = shift @sorted_unique_residue_keys;
        my @rotamer_ids = keys %{ $rotamer_pairs->{$unique_residue_key} };
        my @neighbour_unique_residue_keys =
            map { $_ }
            sort { scalar( keys %{ $rotamer_pairs->{$a} } ) <=>
                   scalar( keys %{ $rotamer_pairs->{$b} } ) }
            keys %{ $residue_pairs->{$unique_residue_key} };

        for my $neighbour_unique_residue_key ( @neighbour_unique_residue_keys ){
            next if $visited_residue_pairs{$unique_residue_key}
                                          {$neighbour_unique_residue_key} ||
                    $visited_residue_pairs{$neighbour_unique_residue_key}
                                          {$unique_residue_key};

            $visited_residue_pairs{$unique_residue_key}
                                  {$neighbour_unique_residue_key} = 1;
            $visited_residue_pairs{$neighbour_unique_residue_key}
                                  {$unique_residue_key} = 1;

            my %pairwise_rotamer_count = ();
            my @neighbour_rotamer_ids =
                keys %{ $rotamer_pairs->{$neighbour_unique_residue_key} };

            for my $rotamer_id ( @rotamer_ids ) {
                if( ! exists $rotamer_atom_site->{$rotamer_id} ) {
                    my %angles =
                        map { $rotamer_angles->{$rotamer_id}{$_}{'type'} =>
                              $rotamer_angles->{$rotamer_id}{$_}{'value'} }
                        keys %{ $rotamer_angles->{$rotamer_id} };
                    my %rotamer_site =
                        %{ clone( $residue_atom_site->{$unique_residue_key} ) };
                    replace_with_rotamer( $parameters, \%rotamer_site,
                                          $unique_residue_key, \%angles );
                    $rotamer_atom_site->{$rotamer_id} = { %rotamer_site };
                }

                for my $neighbour_rotamer_id ( @neighbour_rotamer_ids ) {
                    if( ! exists $rotamer_atom_site->{$neighbour_rotamer_id} &&
                        ! $dry_run ) {
                        my %neighbour_angles =
                            map { $rotamer_angles->{$neighbour_rotamer_id}{$_}{'type'} =>
                                  $rotamer_angles->{$neighbour_rotamer_id}{$_}{'value'} }
                            keys %{ $rotamer_angles->{$neighbour_rotamer_id} };
                        my %neighbour_rotamer_site =
                            %{ clone( $residue_atom_site->{$neighbour_unique_residue_key} ) };
                        replace_with_rotamer( $parameters,
                                              \%neighbour_rotamer_site,
                                              $neighbour_unique_residue_key,
                                              \%neighbour_angles );
                        $rotamer_atom_site->{$neighbour_rotamer_id} =
                            { %neighbour_rotamer_site };
                    }

                    # Calculate pairwise energy.
                    my $pairwise_energy_sum =
                        $dry_run ?
                        0.0 :
                        pairwise_rotamer_energy(
                              $parameters,
                              $rotamer_atom_site->{$rotamer_id},
                              $rotamer_atom_site->{$neighbour_rotamer_id},
                              \&ForceField::NonBonded::general ) / 2;

                    # Under the cutoff limit.
                    if( $pairwise_energy_sum <= $cutoff_atom  ) {
                        $self->{'rotamer_energies'}
                               {$rotamer_id}
                               {$neighbour_rotamer_id} =
                            $pairwise_energy_sum;
                        $self->{'rotamer_energies'}
                               {$neighbour_rotamer_id}
                               {$rotamer_id} =
                            $pairwise_energy_sum;
                        $pairwise_rotamer_count{$rotamer_id} += 1;
                        $pairwise_rotamer_count{$neighbour_rotamer_id} += 1;
                    } else {
                        $pairwise_rotamer_count{$rotamer_id} += 0;
                        $pairwise_rotamer_count{$neighbour_rotamer_id} += 0;
                    }
                }
            }

            # NOTE: need refactoring. Too much loops.
            for my $rotamer_id ( sort keys %pairwise_rotamer_count ) {
                if( ! $pairwise_rotamer_count{$rotamer_id} ) {
                    delete $residue_pairs->{$unique_residue_key}{$rotamer_id};
                }
            }

            if( $verbose ){
                my $rotamer_count =
                    scalar( keys %{ $rotamer_pairs->{$unique_residue_key} } );
                my $neighbour_rotamer_count =
                    scalar( keys %{ $rotamer_pairs->{$neighbour_unique_residue_key} } );

                print info(
                    { message =>
                          "rotamer pairs: " .
                          $unique_residue_key . " " .
                          $rotamer_count . " " .
                          $neighbour_unique_residue_key . " " .
                          $neighbour_rotamer_count . " " .
                          $rotamer_count * $neighbour_rotamer_count . "\n",
                      program => $program_called_by }
                );
            }
        }

        if( $verbose ) {
            # Total count of rotamers for reach residue.
            for my $current_unique_residue_key ( sort keys %{ $rotamer_pairs } ) {
                print info(
                    { message =>
                          "rotamer count: " .
                          $step_count . " " .
                          $current_unique_residue_key . " " .
                          scalar( keys %{ $rotamer_pairs->{$current_unique_residue_key} } ) . "\n",
                      program => $program_called_by }
                )
            }
        }

        $step_count++;
    }

    return;
}

sub get_pairwise_rotamer_energies
{
    my ( $self ) = @_;
    my ( $rotamer_energies ) = ( $self->{'rotamer_energies'} );

    my %pdbx_data = ();
    $pdbx_data{'_[local]_pairwise_energy'}{'metadata'}{'attributes'} =
        [ 'id', 'rotamer_id_1', 'rotamer_id_2', 'type', 'value' ];
    $pdbx_data{'_[local]_pairwise_energy'}{'metadata'}{'is_loop'} = 1;
    $pdbx_data{'_[local]_pairwise_energy'}{'metadata'}{'type'} = 'record';

    my %visited_rotamer_pairs = ();
    my $id = 1;
    for my $rotamer_id (
        sort { $a <=> $b }
        keys %{ $rotamer_energies } ) {
        for my $neighbour_rotamer_id (
            sort { $a <=> $b }
            keys %{ $rotamer_energies->{$rotamer_id} } ) {
            next if $visited_rotamer_pairs{$rotamer_id}{$neighbour_rotamer_id}||
                $visited_rotamer_pairs{$neighbour_rotamer_id}{$rotamer_id};

            $visited_rotamer_pairs{$rotamer_id}{$neighbour_rotamer_id} = 1;
            $visited_rotamer_pairs{$neighbour_rotamer_id}{$rotamer_id} = 1;

            push @{ $pdbx_data{'_[local]_pairwise_energy'}{'data'} },
                { 'id' => $id,
                  'rotamer_id_1' => $rotamer_id,
                  'rotamer_id_2' => $neighbour_rotamer_id,
                  # HACK: for now, the value is hardcoded until the options are
                  # added.
                  'type' => 'composite',
                  'value' =>
                      $rotamer_energies->{$rotamer_id}{$neighbour_rotamer_id} };

            $id++;
        }
    }

    return \%pdbx_data;
}

sub get_rotamer_angles
{
    my ( $self ) = @_;
    my ( $rotamer_angles ) = ( $self->{'rotamer_angles'} );
    my %pdbx_data = ();
    $pdbx_data{'_[local]_rotamer_angle'}{'metadata'}{'attributes'} =
        [ 'id', 'rotamer_id', 'label_seq_id', 'label_comp_id', 'label_asym_id',
          'pdbx_PDB_model_num', 'label_alt_id', 'frequency', 'type', 'value',
          'units' ];
    $pdbx_data{'_[local]_rotamer_angle'}{'metadata'}{'is_loop'} = 1;
    $pdbx_data{'_[local]_rotamer_angle'}{'metadata'}{'type'} = 'record';
    # $pdbx_data{'_[local]_rotamer_angle'}{'data'} = $rotamer_angles;
    return \%pdbx_data;
}

1;
