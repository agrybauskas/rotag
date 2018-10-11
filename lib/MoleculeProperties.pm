package MoleculeProperties;

use strict;
use warnings;

use Exporter qw( import );
our @EXPORT_OK = qw( %RESIDUE_ATOMS );

use Version qw( $VERSION );

our $VERSION = $VERSION;

my %RESIDUE_ATOMS = (
    'SER' => {
        'mandatory' =>
            [
              'N', 'CA', 'C', 'O', 'CB', 'OG',
            ],
        'optional'  =>
            [
              'OXT', 'H', 'H2', 'HA', 'HXT', 'HB2', 'HB3', 'HG',
            ],
    },
    'ARG' => {
        'mandatory' =>
            [
              'N', 'CA', 'C', 'O', 'CB', 'CG', 'CD', 'NE', 'CZ', 'NH1', 'NH2',
            ],
        'optional'  =>
            [
              'OXT', 'H', 'H2', 'HA', 'HXT', 'HB2', 'HB3', 'HG2', 'HG3', 'HD2',
              'HD3', 'HE', 'HH11', 'HH12', 'HH21', 'HH22',
            ],
    },
    'HIS' => {
        'mandatory' =>
            [
              'N', 'CA', 'C', 'O', 'CB', 'CG', 'ND1', 'CD2', 'CE1', 'NE2',
            ],
        'optional'  =>
            [
              'OXT', 'H', 'H2', 'HA', 'HXT', 'HB2', 'HB3', 'HD1', 'HD2', 'HE1',
              'HE2',
            ],
    },
    'LYS' => {
        'mandatory' =>
            [
              'N', 'CA', 'C', 'O', 'CB', 'CG', 'CD', 'CE', 'NZ',
            ],
        'optional'  =>
            [
              'OXT', 'H', 'H2', 'HA', 'HXT', 'HB2', 'HB3', 'HG2', 'HG3', 'HD2',
              'HD3', 'HE2', 'HE3', 'HZ1', 'HZ2',
            ],
    },
    'ASP' => {
        'mandatory' =>
            [
              'N', 'CA', 'C', 'O', 'CB', 'CG', 'OD1', 'OD2',
            ],
        'optional'  =>
            [
              'OXT', 'H', 'H2', 'HA', 'HXT', 'HB2', 'HB3', 'HD2',
            ],
    },
    'GLU' => {
        'mandatory' =>
            [
              'N', 'CA', 'C', 'O', 'CB', 'CG', 'CD', 'OE1', 'OE2',
            ],
        'optional'  =>
            [
              'OXT', 'H', 'H2', 'HA', 'HXT', 'HB2', 'HB3', 'HG2', 'HG3', 'HE2',
            ],
    },
    'CYS' => {
        'mandatory' =>
            [
              'N', 'CA', 'C', 'O', 'CB', 'SG',
            ],
        'optional'  =>
            [
              'OXT', 'H', 'H2', 'HA', 'HXT', 'HB2', 'HB3', 'HG',
            ],
    },
    'GLY' => {
        'mandatory' =>
            [
              'N', 'CA', 'C', 'O',
            ],
        'optional'  =>
            [
              'OXT', 'H', 'H2', 'HA', 'HXT', 'HA2', 'HA3',
            ],
    },
    'PRO' => {
        'mandatory' =>
            [
              'N', 'CA', 'C', 'O', 'CB', 'CG', 'CD',
            ],
        'optional'  =>
            [
              'OXT', 'H', 'HA', 'HXT', 'HB2', 'HB3', 'HG2', 'HG3', 'HD2', 'HD3',
            ],
    },
    'ALA' => {
        'mandatory' =>
            [
              'N', 'CA', 'C', 'O', 'CB',
            ],
        'optional'  =>
            [
              'OXT', 'H', 'H2', 'HA', 'HXT', 'HB1', 'HB2', 'HB3',
            ],
    },
    'VAL' => {
        'mandatory' =>
            [
              'N', 'CA', 'C', 'O', 'CB', 'CG1', 'CG2',
            ],
        'optional'  =>
            [
              'OXT', 'H', 'H2', 'HA', 'HXT', 'HB', 'HG11', 'HG12', 'HG13',
              'HG21', 'HG22', 'HG23',
            ],
    },
    'ILE' => {
        'mandatory' =>
            [
              'N', 'CA', 'C', 'O', 'CB', 'CG1', 'CG2', 'CD1',
            ],
        'optional'  =>
            [
              'OXT', 'H', 'H2', 'HA', 'HXT', 'HB', 'HG12', 'HG13', 'HG21',
              'HG22', 'HG23', 'HD11', 'HD12', 'HD13',
            ],
    },
    'LEU' => {
        'mandatory' =>
            [
              'N', 'CA', 'C', 'O', 'CB', 'CG', 'CD1', 'CD2',
            ],
        'optional'  =>
            [
              'OXT', 'H', 'H2', 'HA', 'HXT','HB2', 'HB3', 'HG', 'HD11', 'HD12',
              'HD13', 'HD21', 'HD22', 'HD23',
            ],
    },
    'MET' => {
        'mandatory' =>
            [
              'N', 'CA', 'C', 'O', 'CB', 'CG', 'SD', 'CE',
            ],
        'optional'  =>
            [
              'OXT', 'H', 'H2', 'HA', 'HXT', 'HB2', 'HB3', 'HG2', 'HG3', 'HE1',
              'HE2', 'HE3',
            ],
    },
    'PHE' => {
        'mandatory' =>
            [
              'N', 'CA', 'C', 'O', 'CB', 'CG', 'CD1', 'CD2', 'CE1', 'CE2', 'CZ',
            ],
        'optional'  =>
            [
              'OXT', 'H', 'H2', 'HA', 'HXT', 'HB2', 'HB3', 'HD1', 'HD2', 'HE1',
              'HE2', 'HZ',
            ],
    },
    'TYR' => {
        'mandatory' =>
            [
              'N', 'CA', 'C', 'O', 'CB', 'CG', 'CD1', 'CD2', 'CE1', 'CE2', 'CZ',
              'OH',
            ],
        'optional'  =>
            [
              'OXT', 'H', 'H2', 'HA', 'HXT', 'HB2', 'HB3', 'HD1', 'HD2', 'HE1',
              'HE2', 'HH',
            ],
    },
    'TRP' => {
        'mandatory' =>
            [
              'N', 'CA', 'C', 'O', 'CB', 'CG', 'CD1', 'CD2', 'NE1', 'CE2',
              'CE3', 'CZ2', 'CZ3', 'CH2',
            ],
        'optional'  =>
            [
              'OXT', 'H', 'H2', 'HA', 'HXT', 'HB2', 'HB3', 'HD1', 'HE1', 'HE3',
              'HZ2', 'HZ3', 'HH2',
            ],
    },
    'THR' => {
        'mandatory' =>
            [
              'N', 'CA', 'C', 'O', 'CB', 'OG1', 'CG2',
            ],
        'optional'  =>
            [
              'OXT', 'H', 'H2', 'HA', 'HXT', 'HB', 'HG1', 'HG21', 'HG22', 'HG23',
            ],
    },
    'ASN' => {
        'mandatory' =>
            [
              'N', 'CA', 'C', 'O', 'CB', 'CG', 'OD1', 'ND2',
            ],
        'optional'  =>
            [
              'OXT', 'H', 'H2', 'HA', 'HXT', 'HB2', 'HB3', 'HD21', 'HD22',
            ],
    },
    'GLN' => {
        'mandatory' =>
            [
              'N', 'CA', 'C', 'O', 'CB', 'CG', 'OD1', 'ND2',
            ],
        'optional'  =>
            [
              'OXT', 'H', 'H2', 'HA', 'HXT', 'HB2', 'HB3', 'HG2', 'HG3',
              'HE21', 'HE22',
            ],
    },
);

1;
