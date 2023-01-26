package BondParameters;

use strict;
use warnings;

use Carp;

use BondPath;
use BondProperties qw( contains_hetatoms
                       contains_sidechain_atoms );
use Measure qw( dihedral_angle );
use PDBxParser qw( expand
                   filter_new
                   split_by );

# ------------------------- Constructors/Destructors -------------------------- #

sub new
{
    my ( $class, $parameters, $atom_site, $options ) = @_;
    my ( $include_mainchain, $include_hetatoms ) =
        ( $options->{'include_mainchain'}, $options->{'include_hetatoms'} );

    $include_mainchain //= 0;
    $include_hetatoms //= 0;

    my $self = { 'parameters' => $parameters,
                 'atom_site'  => $atom_site,
                 'include_mainchain' => $include_mainchain,
                 'include_hetatoms' => $include_hetatoms };

    return bless $self, $class;
}

# ----------------------------- Setters/Getters ------------------------------- #

#
# Identifies bonds that can be rotated by torsional angle.
# Input:
#     $atom_site - atom site data structure (see PDBxParser.pm);
#     $start_atom_ids - starting atom ids;
#     $self->{'include_mainchain'} - flag that includes main-chain atoms;
#     $self->{'include_hetatoms'} - flag that includes heteroatoms.
# Output:
#     data structure that describes rotatable bonds and the constituent atom ids
#     of the bond.
#

sub set_rotatable_bonds
{
    my ( $self, $start_atom_ids ) = @_;
    my ( $include_mainchain, $include_hetatoms ) = (
        $self->{'include_mainchain'},
        $self->{'include_hetatoms'}
    );

    $include_mainchain //= 0;
    $include_hetatoms //= 0;

    my $parameters = $self->{'parameters'};
    my $atom_site = $self->{'atom_site'};

    my $explicit_dihedral_names =
        $self->{'parameters'}{'_[local]_dihedral_angle_name'};

    my $ignore_connections = {
        'label_atom_id' => {
            'N' => { 'CD' => 1, # For PRO.
                     'C' => 1 },
            'C' => { 'O' => 1,
                     ( $include_hetatoms ? ( 'CA' => 1 ) : () ) },
        },
    };

    if( $include_hetatoms ) {
        $start_atom_ids = filter_new(
            $atom_site,
            { 'include' => { 'label_atom_id' => [ 'C' ] },
              'return_data' => 'id'
        } );
    }

    my $bond_paths = BondPath->new( {
        'atom_site' => $atom_site,
        'start_atom_ids' => $start_atom_ids,
        'include_hetatoms' => $include_hetatoms,
        'ignore_connections' => $ignore_connections,
    } );

    my %rotatable_bonds = ();
    my %rotatable_bonds_cache = ();
    my %bond_order = ();
    my $bond_order_idx = 1;

    for my $fourth_atom_id ( @{ $bond_paths->get_atom_order } ) {
        my $third_atom_id = $bond_paths->get_atom_id_to( $fourth_atom_id );

        next if ! defined $third_atom_id;

        my $second_atom_id = $bond_paths->get_atom_id_to( $third_atom_id );

        next if ! defined $second_atom_id;

        my $first_atom_id = $bond_paths->get_atom_id_to( $second_atom_id );

        next if ! defined $first_atom_id;

        # Checks for mainchains and heteroatoms.
        next if ! $include_mainchain &&
            ! contains_sidechain_atoms( $parameters,
                                        $atom_site,
                                        [ $first_atom_id, $second_atom_id,
                                          $third_atom_id, $fourth_atom_id ] ) &&
            ! contains_hetatoms( $atom_site,
                                 [ $first_atom_id, $second_atom_id,
                                   $third_atom_id, $fourth_atom_id ] );

        # Check on hybridization.
        if( ! exists $atom_site->{$second_atom_id}{'hybridization'} ) {
            confess "atom with id $second_atom_id lacks information about " .
                "hybridization";
        }
        if( ! exists $atom_site->{$third_atom_id}{'hybridization'} ){
            confess "atom with id $third_atom_id lacks information " .
                "about hybridization";
        }

        # If one of the bond atoms are sp3, it is rotatable.
        if( $atom_site->{$second_atom_id}{'hybridization'} eq 'sp3' ||
            $atom_site->{$third_atom_id}{'hybridization'} eq 'sp3' ||
            ( $include_hetatoms &&
              $atom_site->{$fourth_atom_id}{'group_PDB'} eq 'HETATM' &&
              $atom_site->{$fourth_atom_id}{'hybridization'} eq '.' ) ) {
            if( exists $rotatable_bonds_cache{$second_atom_id}{$third_atom_id} ) {
                push @{ $rotatable_bonds{$fourth_atom_id} },
                    $rotatable_bonds_cache{$second_atom_id}{$third_atom_id};
            } else {
                push @{ $rotatable_bonds{$fourth_atom_id} },
                    [ $first_atom_id, $second_atom_id, $third_atom_id,
                      $fourth_atom_id ];
                $rotatable_bonds_cache{$second_atom_id}{$third_atom_id} =
                    [ $first_atom_id, $second_atom_id, $third_atom_id,
                      $fourth_atom_id ];
            }
        }

        # If bond atoms are sp2/sp (do not rotate) or just is a continuation of
        # the bond chain, inherits its previous atom's rotatable bonds.
        if( exists $rotatable_bonds{$third_atom_id} ) {
            unshift @{ $rotatable_bonds{$fourth_atom_id} },
                @{ $rotatable_bonds{$third_atom_id} };
        }

        $bond_order{$second_atom_id}{$third_atom_id} = $bond_order_idx;
        $bond_order_idx++;
    }

    # Naming the rotatable bonds.
    my %named_rotatable_bonds = ();
    for my $atom_id ( keys %rotatable_bonds ) {
        my $residue_name = $atom_site->{$atom_id}{'label_comp_id'};
        for my $bond_atom_ids ( @{ $rotatable_bonds{$atom_id} } ) {
            my @rotatable_bond_name_keys =
                ( join( '-', map { $atom_site->{$_}{'label_atom_id'} }
                                @{ $bond_atom_ids } ),
                  join( '-', ( '.',
                               $atom_site->{$bond_atom_ids->[1]}{'label_atom_id'},
                               $atom_site->{$bond_atom_ids->[2]}{'label_atom_id'},
                               '.' ) ) );
            my ( $rotatable_bond_name ) =
                grep { defined $_ }
                     ( ( map { $explicit_dihedral_names->{$residue_name}{$_} }
                             @rotatable_bond_name_keys ),
                       ( map { $explicit_dihedral_names->{'.'}{$_} }
                             @rotatable_bond_name_keys ),
                       $rotatable_bond_name_keys[0] );

            $self->{'dihedral_angles'}{'id'}{$atom_id}{$rotatable_bond_name} = {
                'order' => $bond_order{$bond_atom_ids->[1]}{$bond_atom_ids->[2]},
                'atom_ids' => $bond_atom_ids,
            };
        }
    }

    return;
}

sub set_stretchable_bonds
{

}

sub set_bendable_angles
{

}

sub get_rotatable_bonds
{
    my ( $self ) = @_;
    return $self->{'dihedral_angles'}{'id'};
}

sub get_rotatable_bond_all_atom_ids
{
    my ( $self ) = @_;
    return [ keys %{ $self->{'dihedral_angles'}{'id'} } ];
}

sub get_rotatable_bond_atom_ids
{
    my ( $self, $atom_id, $rotatable_bond_name ) = @_;
    return $self->{'dihedral_angles'}{'id'}{$atom_id}{$rotatable_bond_name}{'atom_ids'};
}

sub get_rotatable_bond_name_order
{
    my ( $self, $atom_id ) = @_;
    return [
        sort { $self->{'dihedral_angles'}{'id'}{$atom_id}{$a}{'order'} <=>
               $self->{'dihedral_angles'}{'id'}{$atom_id}{$b}{'order'} }
        keys %{ $self->{'dihedral_angles'}{'id'}{$atom_id} }
    ];
}

sub get_unique_residue_keys
{
    my ( $self ) = @_;
    return [] if ! exists $self->{'dihedral_angles'}{'unique_residue_key'};
    return [ keys %{ $self->{'dihedral_angles'}{'unique_residue_key'} } ];
}

sub get_dihedral_angles
{
    my ( $self, $unique_residue_key ) = @_;
    return {} if ! exists $self->{'dihedral_angles'}{'unique_residue_key'}
                                                    {$unique_residue_key};
    return $self->{'dihedral_angles'}{'unique_residue_key'}{$unique_residue_key};
}

sub get_dihedral_angle_value
{
    my ( $self, $unique_residue_key, $angle_name ) = @_;
    return $self->{'dihedral_angles'}{'unique_residue_key'}{$unique_residue_key}
                                                           {$angle_name}{'value'};
}

# --------------------------------- Methods ----------------------------------- #

#
# Calculates dihedral angles for all given atoms that are described in atom site
# data structure (produced by obtain_atom_site or functions that uses it). Usage
# of connect_atoms and hybridization functions are necessary for correct
# calculations.
# Input:
#     $atom_site - atom data structure;
#     $self->{'include_mainchain'} - additionally calculates phi and psi
#     mainchain dihedral angles;
#     $self->{'include_hetatoms'} - additionally calculates dihedral angles for
#     hetero atoms.
# Output:
#     data structure that relates residue id and angle values.

sub calculate_dihedral_angles
{
    my ( $self ) = @_;
    my ( $include_mainchain, $include_hetatoms ) = (
        $self->{'include_mainchain'},
        $self->{'include_hetatoms'},
    );

    $include_mainchain //= 0;
    $include_hetatoms //= 0;

    my $parameters = $self->{'parameters'};
    my $atom_site = $self->{'atom_site'};

    my $residue_groups =
        split_by( { 'atom_site' => $atom_site, 'append_dot' => 1 } );

    # Iterates through residue ids and, according to the parameter file,
    # calculates dihedral angles of each rotatable bond.
    for my $residue_unique_key ( sort keys %{ $residue_groups } ) {
        my $residue_site =
            filter_new( $atom_site,
                        { 'include' =>
                          { 'id' => $residue_groups->{$residue_unique_key} } } );

        my $start_atom_ids;
        if( $include_mainchain ) {
            my @expanded_atom_ids = @{ expand( $residue_site, $atom_site, 1 ) };
            my ( $residue_id, $chain_id, $pdbx_model_id, $alt_id ) =
                split ',', $residue_unique_key;
            $start_atom_ids =
                filter_new( $residue_site,
                            { 'include' => { 'id' => \@expanded_atom_ids,
                                             'label_atom_id' => [ 'C' ],
                                             'label_asym_id' => [ $chain_id ],
                                             'pdbx_PDB_model_num' => [ $pdbx_model_id ],
                                             'label_alt_id' => [ $alt_id ] },
                              'exclude' => { 'label_seq_id' => [ $residue_id ] },
                              'return_data' => 'id' } );
            $start_atom_ids = @{ $start_atom_ids } ? $start_atom_ids : undef;
        }

        # NOTE: think about using a some sort of update function.
        if( ! exists $self->{'dihedral_angles'}{'id'} ) {
            $self->set_rotatable_bonds( $start_atom_ids );
        }

        my $unique_rotatable_bonds =
            unique_bond_parameters( $self->{'dihedral_angles'}{'id'} );
        $self->{'dihedral_angles'}{'residue_unique_key'}{$residue_unique_key} =
            $unique_rotatable_bonds;

        # Calculates every side-chain dihedral angle.
        for my $angle_name ( keys %{ $unique_rotatable_bonds } ) {
            my ( $first_atom_id, $second_atom_id, $third_atom_id, $fourth_atom_id ) =
                @{ $unique_rotatable_bonds->{$angle_name}{'atom_ids'} };

            # Extracts coordinates for dihedral angle calculations.
            my ( $first_atom_coord, $second_atom_coord, $third_atom_coord,
                 $fourth_atom_coord ) =
                map { [ $atom_site->{$_}{'Cartn_x'},
                        $atom_site->{$_}{'Cartn_y'},
                        $atom_site->{$_}{'Cartn_z'} ] }
                    ( $first_atom_id, $second_atom_id, $third_atom_id,
                      $fourth_atom_id );
            $self->{'dihedral_angles'}{'residue_unique_key'}
                   {$residue_unique_key}{$angle_name}{'value'} =
                dihedral_angle( [ $first_atom_coord,
                                  $second_atom_coord,
                                  $third_atom_coord,
                                  $fourth_atom_coord ] );
        }
    }

    return;
}

sub calculate_bond_lengths
{

}

sub calculate_bond_angles
{

}

# ----------------------------- Static functions ------------------------------ #

sub unique_bond_parameters
{
    my ( $bond_parameters ) = @_;
    my %unique_bond_parameters;
    for my $atom_id ( sort keys %{ $bond_parameters } ) {
        for my $parameter_name ( keys %{ $bond_parameters->{"$atom_id"} } ) {
            if( ! exists $unique_bond_parameters{"$parameter_name"} ) {
                $unique_bond_parameters{"$parameter_name"} =
                    $bond_parameters->{"$atom_id"}{"$parameter_name"};
            }
        }
    }
    return \%unique_bond_parameters;
}

1;
