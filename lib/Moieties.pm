package Moieties;

use strict;
use warnings;

use Exporter qw( import );
our @EXPORT_OK = qw( %ATOMS
                     %SIDECHAINS );

# ---------------------------------- Moieties --------------------------------- #

our %ATOMS = (
    'H' => {
        1 => {
            'group_PDB' => 'ATOM',
            'id' => 1,
            'type_symbol' => '.',
            'label_atom_id' => undef,
            'label_alt_id' => '.',
            'label_comp_id' => undef,
            'label_asym_id' => undef,
            'label_entity_id' => undef,
            'label_seq_id' => undef,
            'Cartn_x' => -54.128,
            'Cartn_y' => -110.543,
            'Cartn_z' => -6.071,
            'auth_seq_id' => undef,
            'auth_comp_id' => undef,
            'auth_asym_id' => undef,
            'auth_atom_id' => undef,
            'pdbx_PDB_model_num' => undef,
        },
        2 => {
            'group_PDB' => 'ATOM',
            'id' => 2,
            'type_symbol' => 'H',
            'label_atom_id' => undef,
            'label_alt_id' => '.',
            'label_comp_id' => undef,
            'label_asym_id' => undef,
            'label_entity_id' => undef,
            'label_seq_id' => undef,
            'Cartn_x' => -53.464,
            'Cartn_y' => -109.890,
            'Cartn_z' => -5.544,
            'auth_seq_id' => undef,
            'auth_comp_id' => undef,
            'auth_asym_id' => undef,
            'auth_atom_id' => undef,
            'pdbx_PDB_model_num' => undef,
        }
    }
);

our %SIDECHAINS = (
    'SER' => {
        1 => {
            'group_PDB' => 'ATOM',
            'id' => 2,
            'type_symbol' => 'C',
            'label_atom_id' => 'CB',
            'label_alt_id' => '.',
            'label_comp_id' => 'SER',
            'label_asym_id' => undef,
            'label_entity_id' => undef,
            'label_seq_id' => undef,
            'Cartn_x' => -0.494,
            'Cartn_y' => 0.929,
            'Cartn_z' => 0.504,
            'auth_seq_id' => undef,
            'auth_comp_id' => 'SER',
            'auth_asym_id' => undef,
            'auth_atom_id' => 'CB',
            'pdbx_PDB_model_num' => undef,
        },
        2 => {
            'group_PDB' => 'ATOM',
            'id' => 3,
            'type_symbol' => 'O',
            'label_atom_id' => 'OG',
            'label_alt_id' => '.',
            'label_comp_id' => 'SER',
            'label_asym_id' => undef,
            'label_entity_id' => undef,
            'label_seq_id' => undef,
            'Cartn_x' => -0.029,
            'Cartn_y' => 0.446,
            'Cartn_z' => 1.769,
            'auth_seq_id' => undef,
            'auth_comp_id' => 'SER',
            'auth_asym_id' => undef,
            'auth_atom_id' => 'OG',
            'pdbx_PDB_model_num' => undef,
        }
    }
);

1;
