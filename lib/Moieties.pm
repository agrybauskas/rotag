package Moieties;

use strict;
use warnings;

use Exporter qw( import );
our @EXPORT_OK = qw( %ATOMS
                     %SIDECHAINS
                     replace_with_moiety );

use PDBxParser qw( filter_by_unique_residue_key );
use Version qw( $VERSION );

our $VERSION = $VERSION;

# --------------------------------- Moieties ---------------------------------- #

our %ATOMS = (
    'H' => {
        1 => {
            'group_PDB' => 'ATOM',
            'id' => 1,
            'type_symbol' => 'H',
            'label_alt_id' => q{.},
            'Cartn_x' => -53.464,
            'Cartn_y' => -109.890,
            'Cartn_z' => -5.544,
            'pdbx_PDB_model_num' => 1,
        },
    },
);

our %SIDECHAINS = (
    'SER' => {
        1 => {
            'group_PDB' => 'ATOM',
            'id' => 2,
            'type_symbol' => 'C',
            'label_atom_id' => 'CA',
            'label_alt_id' => q{.},
            'label_comp_id' => 'SER',
            'label_asym_id' => 'A',
            'label_entity_id' => 1,
            'label_seq_id' => q{.},
            'Cartn_x' => 0.009,
            'Cartn_y' => 0.077,
            'Cartn_z' => -0.688,
        },
        2 => {
            'group_PDB' => 'ATOM',
            'id' => 2,
            'type_symbol' => 'C',
            'label_atom_id' => 'CB',
            'label_alt_id' => q{.},
            'label_comp_id' => 'SER',
            'label_asym_id' => 'A',
            'label_entity_id' => 1,
            'label_seq_id' => q{.},
            'Cartn_x' => -0.494,
            'Cartn_y' => 0.929,
            'Cartn_z' => 0.504,
        },
        3 => {
            'group_PDB' => 'ATOM',
            'id' => 3,
            'type_symbol' => 'O',
            'label_atom_id' => 'OG',
            'label_alt_id' => q{.},
            'label_comp_id' => 'SER',
            'label_asym_id' => 'A',
            'label_entity_id' => 1,
            'label_seq_id' => q{.},
            'Cartn_x' => -0.029,
            'Cartn_y' => 0.446,
            'Cartn_z' => 1.769,
        },
    },
);

#
# Replaces selected side-chain with specified moiety (usually another
# side-chain).
# Input:
#     $atom_site - atom site data structure (see PDBxParser.pm);
#     $unique_residue_key - key that can determine unique residue
#     (see PDBxParser::unique_residue_key);
#     $moiety - moiety in atom site data structure.
# Output:
#     changes atom site by replacing side-chain with specified moiety.
#

sub replace_with_moiety
{
    my ( $atom_site, $unique_residue_key, $moiety, $options ) = @_;

    my ( $isomer ) = $options->{'isomer'};
    $isomer //= 'R';

    my $residue_site =
        filter_by_unique_residue_key( $atom_site, $unique_residue_key );

}

1;
