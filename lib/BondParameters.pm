package BondParameters;

use strict;
use warnings;

use BondPath;

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
#     $options{'include_mainchain'} - flag that includes main-chain atoms;
#     $options{'include_hetatoms'} - flag that includes heteroatoms.
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

    my $bond_paths = $self->{'bond_paths'};

    $bond_paths //= BondPath->new( {
        'atom_site' => $atom_site,
        'start_atom_ids' => $start_atom_ids,
        'include_hetatoms' => $include_hetatoms,
        'ignore_connections' => $ignore_connections,
    } );
}

#     my %rotatable_bonds = ();
#     my %rotatable_bonds_cache = ();
#     my %bond_order = ();
#     my $bond_order_idx = 1;
#     for my $fourth_atom_id ( sort { $bond_paths->{'atom_order'}{$a} <=>
#                                     $bond_paths->{'atom_order'}{$b} }
#                              keys %{ $bond_paths->{'atom_order'} } ) {

#         my $third_atom_id = $bond_paths->{'to'}{$fourth_atom_id};

#         next if ! defined $third_atom_id;

#         my $second_atom_id = $bond_paths->{'to'}{$third_atom_id};

#         next if ! defined $second_atom_id;

#         my $first_atom_id = $bond_paths->{'to'}{$second_atom_id};

#         next if ! defined $first_atom_id;

#         # Checks for mainchains and heteroatoms.
#         next if ! $include_mainchain &&
#             ! contains_sidechain_atoms( $parameters,
#                                         $atom_site,
#                                         [ $first_atom_id, $second_atom_id,
#                                           $third_atom_id, $fourth_atom_id ] ) &&
#             ! contains_hetatoms( $atom_site,
#                                  [ $first_atom_id, $second_atom_id,
#                                    $third_atom_id, $fourth_atom_id ] );

#         # Check on hybridization.
#         if( ! exists $atom_site->{$second_atom_id}{'hybridization'} ) {
#             confess "atom with id $second_atom_id lacks information about " .
#                 "hybridization";
#         }
#         if( ! exists $atom_site->{$third_atom_id}{'hybridization'} ){
#             confess "atom with id $third_atom_id lacks information " .
#                 "about hybridization";
#         }

#         # If one of the bond atoms are sp3, it is rotatable.
#         if( $atom_site->{$second_atom_id}{'hybridization'} eq 'sp3' ||
#             $atom_site->{$third_atom_id}{'hybridization'} eq 'sp3' ||
#             ( $include_hetatoms &&
#               $atom_site->{$fourth_atom_id}{'group_PDB'} eq 'HETATM' &&
#               $atom_site->{$fourth_atom_id}{'hybridization'} eq '.' ) ) {
#             if( exists $rotatable_bonds_cache{$second_atom_id}{$third_atom_id} ) {
#                 push @{ $rotatable_bonds{$fourth_atom_id} },
#                     $rotatable_bonds_cache{$second_atom_id}{$third_atom_id};
#             } else {
#                 push @{ $rotatable_bonds{$fourth_atom_id} },
#                     [ $first_atom_id, $second_atom_id, $third_atom_id,
#                       $fourth_atom_id ];
#                 $rotatable_bonds_cache{$second_atom_id}{$third_atom_id} =
#                     [ $first_atom_id, $second_atom_id, $third_atom_id,
#                       $fourth_atom_id ];
#             }
#         }

#         # If bond atoms are sp2/sp (do not rotate) or just is a continuation of
#         # the bond chain, inherits its previous atom's rotatable bonds.
#         if( exists $rotatable_bonds{$third_atom_id} ) {
#             unshift @{ $rotatable_bonds{$fourth_atom_id} },
#                 @{ $rotatable_bonds{$third_atom_id} };
#         }

#         $bond_order{$second_atom_id}{$third_atom_id} = $bond_order_idx;
#         $bond_order_idx++;
#     }

#     # Naming the rotatable bonds.
#     my %named_rotatable_bonds = ();
#     for my $atom_id ( keys %rotatable_bonds ) {
#         my $residue_name = $atom_site->{$atom_id}{'label_comp_id'};
#         for my $bond_atom_ids ( @{ $rotatable_bonds{$atom_id} } ) {
#             my @rotatable_bond_name_keys =
#                 ( join( '-', map { $atom_site->{$_}{'label_atom_id'} }
#                                 @{ $bond_atom_ids } ),
#                   join( '-', ( '.',
#                                $atom_site->{$bond_atom_ids->[1]}{'label_atom_id'},
#                                $atom_site->{$bond_atom_ids->[2]}{'label_atom_id'},
#                                '.' ) ) );
#             my ( $rotatable_bond_name ) =
#                 grep { defined $_ }
#                      ( ( map { $explicit_dihedral_names->{$residue_name}{$_} }
#                              @rotatable_bond_name_keys ),
#                        ( map { $explicit_dihedral_names->{'.'}{$_} }
#                              @rotatable_bond_name_keys ),
#                        $rotatable_bond_name_keys[0] );

#             $named_rotatable_bonds{$atom_id}{$rotatable_bond_name} = {
#                 'order' => $bond_order{$bond_atom_ids->[1]}{$bond_atom_ids->[2]},
#                 'atom_ids' => $bond_atom_ids,
#             };
#         }
#     }

#     return \%named_rotatable_bonds;


sub set_stretchable_bonds
{

}

sub set_bendable_angles
{

}

# --------------------------------- Methods ----------------------------------- #

sub calculate_dihedral_angles
{

}

sub calculate_bond_lengths
{

}

sub calculate_bond_angles
{

}

1;
