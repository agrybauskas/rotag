package Moieties;

use strict;
use warnings;

use Exporter qw( import );
our @EXPORT_OK = qw( %SIDECHAINS );

# ---------------------------------- Moieties --------------------------------- #

my %SIDECHAINS = (
    'SER' => {
        'CA' => {
            'group_PDB' = 'ATOM',
            'id' => undef,
            'type_symbol' => 'C',
            'label_atom_id' => 'CA',
            'label_alt_id' => '.',
            'label_comp_id' => 'SER',
            'label_asym_id' => undef,
            'label_entity_id' => undef,
            'label_seq_id' => undef,
            'Cartn_x' => 0.009,
            'Cartn_y' => 0.077,
            'Cartn_z' => -0.688,
            'auth_seq_id' => undef,
            'auth_comp_id' => 'SER',
            'auth_asym_id' => undef,
            'auth_atom_id' => 'CA',
            'pdbx_PDB_model_num' => undef,
        },
        'CB' => {
            'group_PDB' = 'ATOM',
            'id' => undef,
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
        'OG' => {
            'group_PDB' = 'ATOM',
            'id' => undef,
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
