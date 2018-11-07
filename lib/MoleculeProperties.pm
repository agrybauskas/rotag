package MoleculeProperties;

use strict;
use warnings;

use Exporter qw( import );
our @EXPORT_OK = qw( %CLEAR_HYBRIDIZATION
                     %CONNECTIVITY
                     identify_residue_atoms
                     %RESIDUE_ATOMS
                     @ROTATABLE_RESIDUE_NAMES
                     %PARTIAL_CHARGE );

use PDBxParser qw( unique_residue_key
                   split_by );
use Version qw( $VERSION );

our $VERSION = $VERSION;

our %RESIDUE_ATOMS = (
    'XAA' => {
        'mandatory' =>
            [
              'XA', 'XA2', 'CA', 'CB', 'OG',
            ],
        'optional'  =>
            [
              'HA', 'HB2', 'HB3', 'HG', 'CG'
            ],
    },
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
              'OXT', 'H', 'H2', 'HA', 'HXT', 'HB2', 'HB3', 'HD2', 'HE1',
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
              'OXT', 'H', 'H2', 'HXT', 'HA2', 'HA3',
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
              'N', 'CA', 'C', 'O', 'CB', 'CG', 'CD', 'OE1', 'NE2',
            ],
        'optional'  =>
            [
              'OXT', 'H', 'H2', 'HA', 'HXT', 'HB2', 'HB3', 'HG2', 'HG3',
              'HE21', 'HE22',
            ],
    },
);

our %CONNECTIVITY = (
    'XAA' => { # Dummy side-chain.
        'XA'  => [ 'CA' ],
        'XA2' => [ 'CA' ],
        'CA'  => [ 'XA', 'XA2', 'CB' ],
        'CB'  => [ 'CA', 'OG' ],
        'CG'  => [ 'CB' ],
        'OG'  => [ 'CB', 'HA' ],
        'HA'  => [ 'CA' ],
        'HB2' => [ 'CB' ],
        'HB3' => [ 'CB' ],
        'HG'  => [ 'OG' ],
    },
    'SER' => {
        'N'   => [ 'CA', 'H', 'H2' ],
        'CA'  => [ 'N', 'C', 'CB', 'HA' ],
        'C'   => [ 'CA', 'O', 'OXT' ],
        'O'   => [ 'C' ],
        'CB'  => [ 'CA', 'OG', 'HB2', 'HB3' ],
        'OG'  => [ 'CB', 'HG' ],
        'OXT' => [ 'C', 'HXT' ],
        'H'   => [ 'N' ],
        'H2'  => [ 'N' ],
        'HA'  => [ 'CA' ],
        'HXT' => [ 'OXT' ],
        'HB2' => [ 'CB' ],
        'HB3' => [ 'CB' ],
        'HG'  => [ 'OG' ],
    },
    'ARG' => {
        'N'   => [ 'CA', 'H', 'H2' ],
        'CA'  => [ 'N', 'C', 'CB', 'HA' ],
        'C'   => [ 'CA', 'O', 'OXT' ],
        'O'   => [ 'C' ],
        'CB'  => [ 'CA', 'CG', 'HB2', 'HB3' ],
        'CG'  => [ 'CB', 'CD', 'HG2', 'HG3' ],
        'CD'  => [ 'CG', 'NE', 'HD2', 'HD3' ],
        'NE'  => [ 'CD', 'CZ', 'HE' ],
        'CZ'  => [ 'NE', 'NH1', 'NH2' ],
        'NH1' => [ 'CZ', 'HH11', 'HH12' ],
        'NH2' => [ 'CZ', 'HH21', 'HH22' ],
        'OXT' => [ 'C', 'HXT' ],
        'H'   => [ 'N' ],
        'H2'  => [ 'N' ],
        'HA'  => [ 'CA' ],
        'HXT' => [ 'OXT' ],
        'HB2' => [ 'CB' ],
        'HB3' => [ 'CB' ],
        'HG2' => [ 'CG' ],
        'HG3' => [ 'CG' ],
        'HD2' => [ 'CD' ],
        'HD3' => [ 'CD' ],
        'HE'  => [ 'NE' ],
        'HH11'=> [ 'NH1' ],
        'HH12'=> [ 'NH1' ],
        'HH21'=> [ 'NH2' ],
        'HH22'=> [ 'NH2' ],
    },
    'HIS' => {
        'N'   => [ 'CA', 'H', 'H2' ],
        'CA'  => [ 'N', 'C', 'CB', 'HA' ],
        'C'   => [ 'CA', 'O', 'OXT' ],
        'O'   => [ 'C' ],
        'CB'  => [ 'CA', 'CG', 'HB2', 'HB3' ],
        'CG'  => [ 'CB', 'ND1', 'CD2' ],
        'ND1' => [ 'CG', 'CE1', 'HD1' ],
        'CD2' => [ 'CG', 'NE2', 'HD2' ],
        'CE1' => [ 'ND1', 'NE2', 'HE1' ],
        'NE2' => [ 'CD2', 'CE1', 'HE2' ],
        'OXT' => [ 'C', 'HXT' ],
        'H'   => [ 'N' ],
        'H2'  => [ 'N' ],
        'HA'  => [ 'CA' ],
        'HXT' => [ 'OXT' ],
        'HB2' => [ 'CB' ],
        'HB3' => [ 'CB' ],
        'HD1' => [ 'ND1' ],
        'HD2' => [ 'CD2' ],
        'HE1' => [ 'CE1' ],
        'HE2' => [ 'NE2' ],
    },
    'LYS' => {
        'N'   => [ 'CA', 'H', 'H2' ],
        'CA'  => [ 'N', 'C', 'CB', 'HA' ],
        'C'   => [ 'CA', 'O', 'OXT' ],
        'O'   => [ 'C' ],
        'CB'  => [ 'CA', 'CG', 'HB2', 'HB3' ],
        'CG'  => [ 'CB', 'CD', 'HG2', 'HG3' ],
        'CD'  => [ 'CG', 'CE', 'HD2', 'HD3' ],
        'CE'  => [ 'CD', 'NZ', 'HE2', 'HE3' ],
        'NZ'  => [ 'CE', 'HZ1', 'HZ2' ],
        'OXT' => [ 'C', 'HXT' ],
        'H'   => [ 'N' ],
        'H2'  => [ 'N' ],
        'HA'  => [ 'CA' ],
        'HXT' => [ 'OXT' ],
        'HB2' => [ 'CB' ],
        'HB3' => [ 'CB' ],
        'HG2' => [ 'CG' ],
        'HG3' => [ 'CG' ],
        'HD2' => [ 'CD' ],
        'HD3' => [ 'CD' ],
        'HE2' => [ 'CE' ],
        'HE3' => [ 'CE' ],
        'HZ1' => [ 'NZ' ],
        'HZ2' => [ 'NZ' ],
    },
    'ASP' => {
        'N'   => [ 'CA', 'H', 'H2' ],
        'CA'  => [ 'N', 'C', 'CB', 'HA' ],
        'C'   => [ 'CA', 'O', 'OXT' ],
        'O'   => [ 'C' ],
        'CB'  => [ 'CA', 'CG', 'HB2', 'HB3' ],
        'CG'  => [ 'CB', 'OD1', 'OD2' ],
        'OD1' => [ 'CG' ],
        'OD2' => [ 'CG', 'HD2' ],
        'OXT' => [ 'C', 'HXT' ],
        'H'   => [ 'N' ],
        'H2'  => [ 'N' ],
        'HA'  => [ 'CA' ],
        'HXT' => [ 'OXT' ],
        'HB2' => [ 'CB' ],
        'HB3' => [ 'CB' ],
        'HD2' => [ 'OD2' ],
    },
    'GLU' => {
        'N'   => [ 'CA', 'H', 'H2' ],
        'CA'  => [ 'N', 'C', 'CB', 'HA' ],
        'C'   => [ 'CA', 'O', 'OXT' ],
        'O'   => [ 'C' ],
        'CB'  => [ 'CA', 'CG', 'HB2', 'HB3' ],
        'CG'  => [ 'CB', 'CD', 'HG2', 'HG3' ],
        'CD'  => [ 'CG', 'OE1', 'OE2' ],
        'OE1' => [ 'CD' ],
        'OE2' => [ 'CD', 'HE2' ],
        'OXT' => [ 'C', 'HXT' ],
        'H'   => [ 'N' ],
        'H2'  => [ 'N' ],
        'HA'  => [ 'CA' ],
        'HXT' => [ 'OXT' ],
        'HB2' => [ 'CB' ],
        'HB3' => [ 'CB' ],
        'HG2' => [ 'CG' ],
        'HG3' => [ 'CG' ],
        'HE2' => [ 'OE2' ],
    },
    'CYS' => {
        'N'   => [ 'CA', 'H', 'H2' ],
        'CA'  => [ 'N', 'C', 'CB', 'HA' ],
        'C'   => [ 'CA', 'O', 'OXT' ],
        'O'   => [ 'C' ],
        'CB'  => [ 'CA', 'SG', 'HB2', 'HB3' ],
        'SG'  => [ 'CB', 'HG' ],
        'OXT' => [ 'C', 'HXT' ],
        'H'   => [ 'N' ],
        'H2'  => [ 'N' ],
        'HA'  => [ 'CA' ],
        'HXT' => [ 'OXT' ],
        'HB2' => [ 'CB' ],
        'HB3' => [ 'CB' ],
        'HG'  => [ 'SG' ],
    },
    'GLY' => {
        'N'   => [ 'CA', 'H', 'H2' ],
        'CA'  => [ 'N', 'C', 'HA2', 'HA3' ],
        'C'   => [ 'CA', 'O', 'OXT' ],
        'O'   => [ 'C' ],
        'OXT' => [ 'C', 'HXT' ],
        'H'   => [ 'N' ],
        'H2'  => [ 'N' ],
        'HA2'  => [ 'CA' ],
        'HA3'  => [ 'CA' ],
        'HXT' => [ 'OXT' ],
    },
    'PRO' => {
        'N'   => [ 'CA', 'H', 'CD' ],
        'CA'  => [ 'N', 'C', 'CB', 'HA' ],
        'C'   => [ 'CA', 'O', 'OXT' ],
        'O'   => [ 'C' ],
        'CB'  => [ 'CA', 'CG', 'HB2', 'HB3' ],
        'CG'  => [ 'CB', 'CD', 'HG2', 'HG3' ],
        'CD'  => [ 'CG', 'N', 'HD2', 'HD3' ],
        'OXT' => [ 'C', 'HXT' ],
        'H'   => [ 'N' ],
        'HA'  => [ 'CA' ],
        'HXT' => [ 'OXT' ],
        'HB2' => [ 'CB' ],
        'HB3' => [ 'CB' ],
        'HG2' => [ 'CG' ],
        'HG3' => [ 'CG' ],
        'HD2' => [ 'CD' ],
        'HD3' => [ 'CD' ],
    },
    'ALA' => {
        'N'   => [ 'CA', 'H', 'H2' ],
        'CA'  => [ 'N', 'C', 'CB', 'HA' ],
        'C'   => [ 'CA', 'O', 'OXT' ],
        'O'   => [ 'C' ],
        'CB'  => [ 'CA', 'HB1', 'HB2', 'HB3' ],
        'OXT' => [ 'C', 'HXT' ],
        'H'   => [ 'N' ],
        'H2'  => [ 'N' ],
        'HA'  => [ 'CA' ],
        'HXT' => [ 'OXT' ],
        'HB1' => [ 'CB' ],
        'HB2' => [ 'CB' ],
        'HB3' => [ 'CB' ],
    },
    'VAL' => {
        'N'   => [ 'CA', 'H', 'H2' ],
        'CA'  => [ 'N', 'C', 'CB', 'HA' ],
        'C'   => [ 'CA', 'O', 'OXT' ],
        'O'   => [ 'C' ],
        'CB'  => [ 'CA', 'CG1', 'CG2', 'HB' ],
        'CG1' => [ 'CB', 'HG11', 'HG12', 'HG13' ],
        'CG2' => [ 'CB', 'HG21', 'HG22', 'HG23' ],
        'OXT' => [ 'C', 'HXT' ],
        'H'   => [ 'N' ],
        'H2'  => [ 'N' ],
        'HA'  => [ 'CA' ],
        'HXT' => [ 'OXT' ],
        'HB'  => [ 'CB' ],
        'HG11'=> [ 'CG1' ],
        'HG12'=> [ 'CG1' ],
        'HG13'=> [ 'CG1' ],
        'HG21'=> [ 'CG2' ],
        'HG22'=> [ 'CG2' ],
        'HG23'=> [ 'CG2' ],
    },
    'ILE' => {
        'N'   => [ 'CA', 'H', 'H2' ],
        'CA'  => [ 'N', 'C', 'CB', 'HA' ],
        'C'   => [ 'CA', 'O', 'OXT' ],
        'O'   => [ 'C' ],
        'CB'  => [ 'CA', 'CG1', 'CG2', 'HB' ],
        'CG1' => [ 'CB', 'CD1', 'HG12', 'HG13' ],
        'CG2' => [ 'CB', 'HG21', 'HG22', 'HG23' ],
        'CD1' => [ 'CB', 'HD11', 'HD12', 'HD13' ],
        'OXT' => [ 'C', 'HXT' ],
        'H'   => [ 'N' ],
        'H2'  => [ 'N' ],
        'HA'  => [ 'CA' ],
        'HXT' => [ 'OXT' ],
        'HB'  => [ 'CB' ],
        'HG12'=> [ 'CG1' ],
        'HG13'=> [ 'CG1' ],
        'HG21'=> [ 'CG2' ],
        'HG22'=> [ 'CG2' ],
        'HG23'=> [ 'CG2' ],
        'HD11'=> [ 'CD1' ],
        'HD12'=> [ 'CD1' ],
        'HD13'=> [ 'CD1' ],
    },
    'LEU' => {
        'N'   => [ 'CA', 'H', 'H2' ],
        'CA'  => [ 'N', 'C', 'CB', 'HA' ],
        'C'   => [ 'CA', 'O', 'OXT' ],
        'O'   => [ 'C' ],
        'CB'  => [ 'CA', 'CG', 'HB2', 'HB3' ],
        'CG'  => [ 'CB', 'CD1', 'CD2', 'HG' ],
        'CD1' => [ 'CG', 'HD11', 'HD12', 'HD13' ],
        'CD2' => [ 'CG', 'HD21', 'HD22', 'HD23' ],
        'OXT' => [ 'C', 'HXT' ],
        'H'   => [ 'N' ],
        'H2'  => [ 'N' ],
        'HA'  => [ 'CA' ],
        'HXT' => [ 'OXT' ],
        'HB2' => [ 'CB' ],
        'HB3' => [ 'CB' ],
        'HG'  => [ 'CG' ],
        'HD11'=> [ 'CD1' ],
        'HD12'=> [ 'CD1' ],
        'HD13'=> [ 'CD1' ],
        'HD21'=> [ 'CD2' ],
        'HD22'=> [ 'CD2' ],
        'HD23'=> [ 'CD2' ],
    },
    'MET' => {
        'N'   => [ 'CA', 'H', 'H2' ],
        'CA'  => [ 'N', 'C', 'CB', 'HA' ],
        'C'   => [ 'CA', 'O', 'OXT' ],
        'O'   => [ 'C' ],
        'CB'  => [ 'CA', 'CG', 'HB2', 'HB3' ],
        'CG'  => [ 'CB', 'SD', 'HG2', 'HG3' ],
        'SD'  => [ 'CG', 'CE' ],
        'CE'  => [ 'SD', 'HE1', 'HE2', 'HE3' ],
        'OXT' => [ 'C', 'HXT' ],
        'H'   => [ 'N' ],
        'H2'  => [ 'N' ],
        'HA'  => [ 'CA' ],
        'HXT' => [ 'OXT' ],
        'HB2' => [ 'CB' ],
        'HB3' => [ 'CB' ],
        'HG2' => [ 'CG' ],
        'HG3' => [ 'CG' ],
        'HE1' => [ 'CE' ],
        'HE2' => [ 'CE' ],
        'HE3' => [ 'CE' ],
    },
    'PHE' => {
        'N'   => [ 'CA', 'H', 'H2' ],
        'CA'  => [ 'N', 'C', 'CB', 'HA' ],
        'C'   => [ 'CA', 'O', 'OXT' ],
        'O'   => [ 'C' ],
        'CB'  => [ 'CA', 'CG', 'HB2', 'HB3' ],
        'CG'  => [ 'CB', 'CD1', 'CD2' ],
        'CD1' => [ 'CG', 'CE1', 'HD1' ],
        'CD2' => [ 'CG', 'CE2', 'HD2' ],
        'CE1' => [ 'CD1', 'CZ', 'HE1' ],
        'CE2' => [ 'CD2', 'CZ', 'HE2' ],
        'CZ'  => [ 'CE1', 'CE2', 'HZ' ],
        'OXT' => [ 'C', 'HXT' ],
        'H'   => [ 'N' ],
        'H2'  => [ 'N' ],
        'HA'  => [ 'CA' ],
        'HXT' => [ 'OXT' ],
        'HB2' => [ 'CB' ],
        'HB3' => [ 'CB' ],
        'HD1' => [ 'CD1' ],
        'HD2' => [ 'CD2' ],
        'HE1' => [ 'CE1' ],
        'HE2' => [ 'CE2' ],
        'HZ'  => [ 'CZ' ],
    },
    'TYR' => {
        'N'   => [ 'CA', 'H', 'H2' ],
        'CA'  => [ 'N', 'C', 'CB', 'HA' ],
        'C'   => [ 'CA', 'O', 'OXT' ],
        'O'   => [ 'C' ],
        'CB'  => [ 'CA', 'CG', 'HB2', 'HB3' ],
        'CG'  => [ 'CB', 'CD1', 'CD2' ],
        'CD1' => [ 'CG', 'CE1', 'HD1' ],
        'CD2' => [ 'CG', 'CE2', 'HD2' ],
        'CE1' => [ 'CD1', 'CZ', 'HE1' ],
        'CE2' => [ 'CD2', 'CZ', 'HE2' ],
        'CZ'  => [ 'CE1', 'CE2', 'OH' ],
        'OH'  => [ 'CZ', 'HH' ],
        'OXT' => [ 'C', 'HXT' ],
        'H'   => [ 'N' ],
        'H2'  => [ 'N' ],
        'HA'  => [ 'CA' ],
        'HXT' => [ 'OXT' ],
        'HB2' => [ 'CB' ],
        'HB3' => [ 'CB' ],
        'HD1' => [ 'CD1' ],
        'HD2' => [ 'CD2' ],
        'HE1' => [ 'CE1' ],
        'HE2' => [ 'CE2' ],
        'HH'  => [ 'CZ' ],
    },
    'TRP' => {
        'N'   => [ 'CA', 'H', 'H2' ],
        'CA'  => [ 'N', 'C', 'CB', 'HA' ],
        'C'   => [ 'CA', 'O', 'OXT' ],
        'O'   => [ 'C' ],
        'CB'  => [ 'CA', 'CG', 'HB2', 'HB3' ],
        'CG'  => [ 'CB', 'CD1', 'CD2' ],
        'CD1' => [ 'CG', 'NE1', 'HD1' ],
        'CD2' => [ 'CG', 'CE2', 'CE3' ],
        'NE1' => [ 'CD1', 'CE2', 'HE1' ],
        'CE2' => [ 'CD2', 'NE1', 'CZ2' ],
        'CE3' => [ 'CD2', 'CZ3', 'HE3' ],
        'CZ2' => [ 'CE2', 'CH2', 'HZ2' ],
        'CZ3' => [ 'CE3', 'CH2', 'HZ3' ],
        'CH2' => [ 'CZ2', 'CZ3', 'HH2' ],
        'OXT' => [ 'C', 'HXT' ],
        'H'   => [ 'N' ],
        'H2'  => [ 'N' ],
        'HA'  => [ 'CA' ],
        'HXT' => [ 'OXT' ],
        'HB2' => [ 'CB' ],
        'HB3' => [ 'CB' ],
        'HD1' => [ 'CD1' ],
        'HE1' => [ 'NE1' ],
        'HE3' => [ 'CE3' ],
        'HZ2' => [ 'CZ2' ],
        'HZ3' => [ 'CZ3' ],
        'HH2' => [ 'CH2' ],
    },
    'THR' => {
        'N'   => [ 'CA', 'H', 'H2' ],
        'CA'  => [ 'N', 'C', 'CB', 'HA' ],
        'C'   => [ 'CA', 'O', 'OXT' ],
        'O'   => [ 'C' ],
        'CB'  => [ 'CA', 'OG1', 'CG2', 'HB' ],
        'OG1' => [ 'CB', 'HG1' ],
        'CG2' => [ 'CB', 'HG21', 'HG22', 'HG23' ],
        'OXT' => [ 'C', 'HXT' ],
        'H'   => [ 'N' ],
        'H2'  => [ 'N' ],
        'HA'  => [ 'CA' ],
        'HXT' => [ 'OXT' ],
        'HB'  => [ 'CB' ],
        'HG1'=> [ 'OG1' ],
        'HG21'=> [ 'CG2' ],
        'HG22'=> [ 'CG2' ],
        'HG23'=> [ 'CG2' ],
    },
    'ASN' => {
        'N'   => [ 'CA', 'H', 'H2' ],
        'CA'  => [ 'N', 'C', 'CB', 'HA' ],
        'C'   => [ 'CA', 'O', 'OXT' ],
        'O'   => [ 'C' ],
        'CB'  => [ 'CA', 'CG', 'HB2', 'HB3' ],
        'CG'  => [ 'CB', 'OD1', 'ND2' ],
        'OD1' => [ 'CG' ],
        'ND2' => [ 'CG', 'HD21', 'HD22' ],
        'OXT' => [ 'C', 'HXT' ],
        'H'   => [ 'N' ],
        'H2'  => [ 'N' ],
        'HA'  => [ 'CA' ],
        'HXT' => [ 'OXT' ],
        'HB2' => [ 'CB' ],
        'HB3' => [ 'CB' ],
        'HD21' => [ 'ND2' ],
        'HD22' => [ 'ND2' ],
    },
    'GLN' => {
        'N'   => [ 'CA', 'H', 'H2' ],
        'CA'  => [ 'N', 'C', 'CB', 'HA' ],
        'C'   => [ 'CA', 'O', 'OXT' ],
        'O'   => [ 'C' ],
        'CB'  => [ 'CA', 'CG', 'HB2', 'HB3' ],
        'CG'  => [ 'CB', 'CD', 'HG2', 'HG3' ],
        'CD'  => [ 'CG', 'OE1', 'NE2' ],
        'OE1' => [ 'CD' ],
        'NE2' => [ 'CD', 'HE21', 'HE22' ],
        'OXT' => [ 'C', 'HXT' ],
        'H'   => [ 'N' ],
        'H2'  => [ 'N' ],
        'HA'  => [ 'CA' ],
        'HXT' => [ 'OXT' ],
        'HB2' => [ 'CB' ],
        'HB3' => [ 'CB' ],
        'HG2' => [ 'CG' ],
        'HG3' => [ 'CG' ],
        'HE21' => [ 'NE2' ],
        'HE22' => [ 'NE2' ],
    },
);

# Parameters taken from:
# Duan, Yong, Chun Wu, Shibasish Chowdhury, Mathew C. Lee, Guoming
# Xiong, Wei Zhang, Rong Yang et al. "A point‐charge force field for
# molecular mechanics simulations of proteins based on condensed‐phase
# quantum mechanical calculations." Journal of computational chemistry
# 24, no. 16 (2003): 1999-2012.
#
# TODO: because it is Amber force field, force field parameters should be taken
# from original Amber files.

our %PARTIAL_CHARGE = (
    'XAA' => { # Dummy side-chain.
        'XA'  => 0.0,
        'XA2' => 0.0,
        'CA'  => 0.118,
        'CB'  => 0.147,
        'CG'  => 0.147,
        'OG'  => -0.640,
        'HA'  => 0.142,
        'HB2' => 0.04,
        'HB3' => 0.04,
        'HG'  => 0.446,
    },
    'SER' => {
        'N'   => -0.541,
        'CA'  => 0.118,
        'C'   => 0.483,
        'O'   => -0.581,
        'CB'  => 0.147,
        'OG'  => -0.640,
        'OXT' => -0.581, # HACK: make a correction of true OXT (now same as O).
        'H2'  => 0.345,  # HACK: make a correction of true H2 (now same as H).
        'H'   => 0.345,
        'HA'  => 0.142,
        'HB2' => 0.04,
        'HB3' => 0.04,
        'HG'  => 0.446,
    },
    'ARG' => {
        'N'   => -0.301,
        'CA'  => -0.131,
        'C'   => 0.73,
        'O'   => -0.578,
        'CB'  => 0.037,
        'CG'  => 0.012,
        'CD'  => 0.126,
        'NE'  => 0.465,
        'CZ'  => 0.566,
        'NH1' => -0.686,
        'NH2' => -0.686,
        'OXT' => -0.581, # HACK: make a correction of true OXT (now same as O).
        'H2'  => 0.345,  # HACK: make a correction of true H2 (now same as H).
        'H'   => 0.234,
        'HA'  => 0.053,
        'HB2' => 0.028,
        'HB3' => 0.028,
        'HG2' => 0.003,
        'HG3' => 0.003,
        'HD2' => 0.068,
        'HD3' => 0.068,
        'HE'  => 0.326,
        'HH11'=> 0.391,
        'HH12'=> 0.391,
        'HH21'=> 0.391,
        'HH22'=> 0.391,
    },
    'HIS' => {
        'N'   => -0.528,
        'CA'  => 0.031,
        'C'   => 0.662,
        'O'   => -0.529,
        'CB'  => -0.152,
        'CG'  => 0.278,
        'ND1' => -0.423,
        'CD2' => -0.298,
        'CE1' => 0.026,
        'NE2' => -0.098,
        'OXT' => -0.581, # HACK: make a correction of true OXT (now same as O).
        'H2'  => 0.345,  # HACK: make a correction of true H2 (now same as H).
        'H'   => 0.282,
        'HA'  => 0.085,
        'HB2' => 0.055,
        'HB3' => 0.055,
        'HD2' => 0.16,
        'HE1' => 0.127,
        'HE2' => 0.267,
    },
    'LYS' => {
        'N'   => -0.436,
        'CA'  => -0.039,
        'C'   => 0.725,
        'O'   => -0.563,
        'CB'  => -0.108,
        'CG'  => 0.033,
        'CD'  => -0.048,
        'CE'  => -0.070,
        'NZ'  => -0.250,
        'OXT' => -0.581, # HACK: make a correction of true OXT (now same as O).
        'H2'  => 0.345,  # HACK: make a correction of true H2 (now same as H).
        'H'   => 0.251,
        'HA'  => 0.129,
        'HB2' => 0.045,
        'HB3' => 0.045,
        'HG2' => 0.01,
        'HG3' => 0.01,
        'HD2' => 0.071,
        'HD3' => 0.071,
        'HE2' => 0.12,
        'HE3' => 0.12,
        'HZ1' => 0.295,
        'HZ2' => 0.295,
    },
    'ASP' => {
        'N'   => -0.558,
        'CA'  => 0.007,
        'C'   => 0.443,
        'O'   => -0.501,
        'CB'  => -0.048,
        'CG'  => 0.745,
        'OD1' => -0.730,
        'OD2' => -0.730,
        'OXT' => -0.581, # HACK: make a correction of true OXT (now same as O).
        'H2'  => 0.345,  # HACK: make a correction of true H2 (now same as H).
        'H'   => 0.32,
        'HA'  => 0.082,
        'HB2' => -0.015,
        'HB3' => -0.015,
    },
    'GLU' => {
        'N'   => -0.423,
        'CA'  => 0.032,
        'C'   => 0.47,
        'O'   => -0.593,
        'CB'  => 0.075,
        'CG'  => -0.034,
        'CD'  => 0.765,
        'OE1' => -0.824,
        'OE2' => -0.824,
        'OXT' => -0.581, # HACK: make a correction of true OXT (now same as O).
        'H2'  => 0.345,  # HACK: make a correction of true H2 (now same as H).
        'H'   => 0.307,
        'HA'  => 0.065,
        'HB2' => -0.004,
        'HB3' => -0.004,
        'HG2' => -0.004,
        'HG3' => -0.004,
    },
    'CYS' => {
        'N'   => -0.396,
        'CA'  => -0.074,
        'C'   => 0.643,
        'O'   => -0.585,
        'CB'  => -0.221,
        'SG'  => -0.285,
        'OXT' => -0.581, # HACK: make a correction of true OXT (now same as O).
        'H2'  => 0.345,  # HACK: make a correction of true H2 (now same as H).
        'H'   => 0.295,
        'HA'  => 0.141,
        'HB2' => 0.147,
        'HB3' => 0.147,
        'HG'  => 0.189,
    },
    'GLY' => {
        'N'   => -0.374,
        'CA'  => -0.129,
        'C'   => 0.581,
        'O'   => -0.509,
        'OXT' => -0.581, # HACK: make a correction of true OXT (now same as O).
        'H2'  => 0.345,  # HACK: make a correction of true H2 (now same as H).
        'H'   => 0.254,
        'HA2'  => 0.089,
        'HA3'  => 0.089,
    },
    'PRO' => {
        'N'   => -0.088,
        'CA'  => -0.035,
        'C'   => 0.334,
        'O'   => -0.435,
        'CB'  => -0.003,
        'CG'  => 0.013,
        'CD'  => -0.012,
        'OXT' => -0.581, # HACK: make a correction of true OXT (now same as O).
        'H'   => 0.345,  # HACK: make a correction of true H2 (now same as H).
        'HA'  => 0.06,
        'HB2' => 0.019,
        'HB3' => 0.019,
        'HG2' => 0.02,
        'HG3' => 0.02,
        'HD2' => 0.044,
        'HD3' => 0.044,
    },
    'ALA' => {
        'N'   => -0.405,
        'CA'  => -0.028,
        'C'   => 0.57,
        'O'   => -0.555,
        'CB'  => -0.230,
        'OXT' => -0.581, # HACK: make a correction of true OXT (now same as O).
        'H2'  => 0.345,  # HACK: make a correction of true H2 (now same as H).
        'H'   => 0.294,
        'HA'  => 0.121,
        'HB1' => 0.077,
        'HB2' => 0.077,
        'HB3' => 0.077,
    },
    'VAL' => {
        'N'   => -0.450,
        'CA'  => -0.052,
        'C'   => 0.447,
        'O'   => -0.405,
        'CB'  => 0.395,
        'CG1' => -0.090,
        'CG2' => -0.090,
        'OXT' => -0.581, # HACK: make a correction of true OXT (now same as O).
        'H2'  => 0.345,  # HACK: make a correction of true H2 (now same as H).
        'H'   => 0.44,
        'HA'  => -0.026,
        'HB'  => -0.116,
        'HG11'=> -0.009,
        'HG12'=> -0.009,
        'HG13'=> -0.009,
        'HG21'=> -0.009,
        'HG22'=> -0.009,
        'HG23'=> -0.009,
    },
    'ILE' => {
        'N'   => -0.451,
        'CA'  => -0.102,
        'C'   => 0.569,
        'O'   => -0.620,
        'CB'  => 0.062,
        'CG1' => 0.022,
        'CG2' => -0.130,
        'CD1' => -0.101,
        'OXT' => -0.581, # HACK: make a correction of true OXT (now same as O).
        'H2'  => 0.345,  # HACK: make a correction of true H2 (now same as H).
        'H'   => 0.329,
        'HA'  => 0.174,
        'HB'  => 0.062,
        'HG12'=> 0.012,
        'HG13'=> 0.012,
        'HG21'=> 0.03,
        'HG22'=> 0.03,
        'HG23'=> 0.03,
        'HD11'=> 0.024,
        'HD12'=> 0.024,
        'HD13'=> 0.024,
    },
    'LEU' => {
        'N'   => -0.355,
        'CA'  => -0.101,
        'C'   => 0.573,
        'O'   => -0.558,
        'CB'  => -0.144,
        'CG'  => 0.192,
        'CD1' => -0.123,
        'CD2' => -0.123,
        'OXT' => -0.581, # HACK: make a correction of true OXT (now same as O).
        'H2'  => 0.345,  # HACK: make a correction of true H2 (now same as H).
        'H'   => 0.262,
        'HA'  => 0.137,
        'HB2' => 0.053,
        'HB3' => 0.053,
        'HG'  => 0.001,
        'HD11'=> 0.022,
        'HD12'=> 0.022,
        'HD13'=> 0.022,
        'HD21'=> 0.022,
        'HD22'=> 0.022,
        'HD23'=> 0.022,
    },
    'MET' => {
        'N'   => -0.395,
        'CA'  => -0.088,
        'C'   => 0.6,
        'O'   => -0.566,
        'CB'  => 0.019,
        'CG'  => -0.208,
        'SD'  => -0.212,
        'CE'  => -0.285,
        'OXT' => -0.581, # HACK: make a correction of true OXT (now same as O).
        'H2'  => 0.345,  # HACK: make a correction of true H2 (now same as H).
        'H'   => 0.281,
        'HA'  => 0.123,
        'HB2' => 0.049,
        'HB3' => 0.049,
        'HG2' => 0.124,
        'HG3' => 0.124,
        'HE1' => 0.128,
        'HE2' => 0.128,
        'HE3' => 0.128,
    },
    'PHE' => {
        'N'   => -0.371,
        'CA'  => -0.030,
        'C'   => 0.548,
        'O'   => -0.507,
        'CB'  => -0.099,
        'CG'  => 0.021,
        'CD1' => -0.083,
        'CD2' => -0.083, # NOTE: there was no value in the article.
        'CE1' => -0.157,
        'CE2' => -0.157, # NOTE: there was no value in the article.
        'CZ'  => -0.100,
        'OXT' => -0.581, # HACK: make a correction of true OXT (now same as O).
        'H2'  => 0.345,  # HACK: make a correction of true H2 (now same as H).
        'H'   => 0.234,
        'HA'  => 0.102,
        'HB2' => 0.061,
        'HB3' => 0.061,
        'HD1' => 0.098,
        'HD2' => 0.098,
        'HE1' => 0.124,
        'HE2' => 0.124,
        'HZ'  => 0.115,
    },
    'TYR' => {
        'N'   => -0.488,
        'CA'  => 0.010,
        'C'   => 0.622,
        'O'   => -0.527,
        'CB'  => -0.052,
        'CG'  => 0.113,
        'CD1' => -0.183,
        'CE1' => -0.182,
        'CE2' => -0.182, # NOTE: there was no value in the article.
        'CZ'  => 0.206,
        'OH'  => -0.421,
        'OXT' => -0.581, # HACK: make a correction of true OXT (now same as O).
        'H2'  => 0.345,  # HACK: make a correction of true H2 (now same as H).
        'H'   => 0.264,
        'HA'  => 0.096,
        'HB2' => 0.019,
        'HB3' => 0.019,
        'HD1' => 0.133,
        'HD2' => 0.133,
        'HE1' => 0.137,
        'HE2' => 0.137,
        'HH'  => 0.330,
    },
    'TRP' => {
        'N'   => -0.428,
        'CA'  => -0.020,
        'C'   => 0.584,
        'O'   => -0.495,
        'CB'  => -0.098,
        'CG'  => -0.100,
        'CD1' => -0.174,
        'CD2' => 0.09,
        'NE1' => -0.298,
        'CE2' => 0.142,
        'CE3' => -0.154,
        'CZ2' => -0.211,
        'CZ3' => -0.164,
        'CH2' => -0.133,
        'OXT' => -0.581, # HACK: make a correction of true OXT (now same as O).
        'H2'  => 0.345,  # HACK: make a correction of true H2 (now same as H).
        'H'   => 0.242,
        'HA'  => 0.107,
        'HB2' => 0.065,
        'HB3' => 0.065,
        'HD1' => 0.171,
        'HE1' => 0.322,
        'HE3' => 0.123,
        'HZ2' => 0.126,
        'HZ3' => 0.119,
        'HH2' => 0.119,
    },
    'THR' => {
        'N'   => -0.245,
        'CA'  => -0.271,
        'C'   => 0.56,
        'O'   => -0.552,
        'CB'  => 0.238,
        'OG1' => -0.602,
        'CG2' => -0.176,
        'OXT' => -0.581, # HACK: make a correction of true OXT (now same as O).
        'H2'  => 0.345,  # HACK: make a correction of true H2 (now same as H).
        'H'   => 0.255,
        'HA'  => 0.164,
        'HB'  => 0.045,
        'HG1'=> 0.405,
        'HG21'=> 0.06,
        'HG22'=> 0.06,
        'HG23'=> 0.06,
    },
    'ASN' => {
        'N'   => -0.430,
        'CA'  => 0.045,
        'C'   => 0.617,
        'O'   => -0.524,
        'CB'  => -0.094,
        'CG'  => 0.584,
        'OD1' => -0.527,
        'ND2' => -0.782,
        'OXT' => -0.581, # HACK: make a correction of true OXT (now same as O).
        'H2'  => 0.345,  # HACK: make a correction of true H2 (now same as H).
        'H'   => 0.255,
        'HA'  => 0.06,
        'HB2' => 0.043,
        'HB3' => 0.043,
        'HD21' => 0.355,
        'HD22' => 0.355,
    },
    'GLN' => {
        'N'   => -0.387,
        'CA'  => 0.037,
        'C'   => 0.419,
        'O'   => -0.565,
        'CB'  => -0.032,
        'CG'  => -0.020,
        'CD'  => 0.668,
        'OE1' => -0.628,
        'NE2' => -0.883,
        'OXT' => -0.581, # HACK: make a correction of true OXT (now same as O).
        'H2'  => 0.345,  # HACK: make a correction of true H2 (now same as H).
        'H'   => 0.301,
        'HA'  => 0.152,
        'HB2' => 0.031,
        'HB3' => 0.031,
        'HG2' => -0.031,
        'HG3' => 0.031,
        'HE21' => 0.408,
        'HE22' => 0.408,
    },
);

our %CLEAR_HYBRIDIZATION = (
    'XAA' => { # Dummy side-chain.
        'XA'  => 'sp3',
        'XA2' => 'sp3',
        'CA'  => 'sp3',
        'CB'  => 'sp3',
        'CG'  => 'sp3',
        'OG'  => 'sp3',
    },
    'SER' => {
        'CA'  => 'sp3',
        'C'   => 'sp2',
        'O'   => 'sp2',
        'CB'  => 'sp3',
        'OG'  => 'sp3',
        'OXT' => 'sp2',
    },
    'ARG' => {
        'CA'  => 'sp3',
        'C'   => 'sp2',
        'O'   => 'sp2',
        'CB'  => 'sp3',
        'CG'  => 'sp3',
        'CD'  => 'sp3',
        'NE'  => 'sp2',
        'CZ'  => 'sp2',
        'NH1' => 'sp2',
        'NH2' => 'sp2',
        'OXT' => 'sp2',
    },
    'HIS' => {
        'CA'  => 'sp3',
        'C'   => 'sp2',
        'O'   => 'sp2',
        'CB'  => 'sp3',
        'CG'  => 'sp2',
        'ND1' => 'sp2',
        'CD2' => 'sp2',
        'CE1' => 'sp2',
        'NE2' => 'sp2',
        'OXT' => 'sp2',
    },
    'LYS' => {
        'CA'  => 'sp3',
        'C'   => 'sp2',
        'O'   => 'sp2',
        'CB'  => 'sp3',
        'CG'  => 'sp3',
        'CD'  => 'sp3',
        'CE'  => 'sp3',
        'NZ'  => 'sp3',
        'OXT' => 'sp2',
    },
    'ASP' => {
        'CA'  => 'sp3',
        'C'   => 'sp2',
        'O'   => 'sp2',
        'CB'  => 'sp3',
        'CG'  => 'sp2',
        'OD1' => 'sp2',
        'OD2' => 'sp2',
        'OXT' => 'sp2',
    },
    'GLU' => {
        'CA'  => 'sp3',
        'C'   => 'sp2',
        'O'   => 'sp2',
        'CB'  => 'sp3',
        'CG'  => 'sp3',
        'CD'  => 'sp2',
        'OE1' => 'sp2',
        'OE2' => 'sp2',
        'OXT' => 'sp2',
    },
    'CYS' => {
        'CA'  => 'sp3',
        'C'   => 'sp2',
        'O'   => 'sp2',
        'CB'  => 'sp3',
        'SG'  => 'sp3',
        'OXT' => 'sp2',
    },
    'GLY' => {
        'CA'  => 'sp3',
        'C'   => 'sp2',
        'O'   => 'sp2',
        'OXT' => 'sp2',
    },
    'PRO' => {
        'CA'  => 'sp3',
        'C'   => 'sp2',
        'O'   => 'sp2',
        'CB'  => 'sp3',
        'CG'  => 'sp3',
        'CD'  => 'sp3',
        'OXT' => 'sp2',
    },
    'ALA' => {
        'CA'  => 'sp3',
        'C'   => 'sp2',
        'O'   => 'sp2',
        'CB'  => 'sp3',
        'OXT' => 'sp2',
    },
    'VAL' => {
        'CA'  => 'sp3',
        'C'   => 'sp2',
        'O'   => 'sp2',
        'CB'  => 'sp3',
        'CG1' => 'sp3',
        'CG2' => 'sp3',
        'OXT' => 'sp2',
    },
    'ILE' => {
        'CA'  => 'sp3',
        'C'   => 'sp2',
        'O'   => 'sp2',
        'CB'  => 'sp3',
        'CG1' => 'sp3',
        'CG2' => 'sp3',
        'CD1' => 'sp3',
        'OXT' => 'sp2',
    },
    'LEU' => {
        'CA'  => 'sp3',
        'C'   => 'sp2',
        'O'   => 'sp2',
        'CB'  => 'sp3',
        'CG'  => 'sp3',
        'CD1' => 'sp3',
        'CD2' => 'sp3',
        'OXT' => 'sp2',
    },
    'MET' => {
        'CA'  => 'sp3',
        'C'   => 'sp2',
        'O'   => 'sp2',
        'CB'  => 'sp3',
        'CG'  => 'sp3',
        'SD'  => 'sp3',
        'CE'  => 'sp3',
        'OXT' => 'sp2',
    },
    'PHE' => {
        'CA'  => 'sp3',
        'C'   => 'sp2',
        'O'   => 'sp2',
        'CB'  => 'sp3',
        'CG'  => 'sp2',
        'CD1' => 'sp2',
        'CD2' => 'sp2',
        'CE1' => 'sp2',
        'CE2' => 'sp2',
        'CZ'  => 'sp2',
        'OXT' => 'sp2',
    },
    'TYR' => {
        'CA'  => 'sp3',
        'C'   => 'sp2',
        'O'   => 'sp2',
        'CB'  => 'sp3',
        'CG'  => 'sp2',
        'CD1' => 'sp2',
        'CE1' => 'sp2',
        'CE2' => 'sp2',
        'CZ'  => 'sp2',
        'OH'  => 'sp3',
        'OXT' => 'sp3',
    },
    'TRP' => {
        'CA'  => 'sp3',
        'C'   => 'sp2',
        'O'   => 'sp2',
        'CB'  => 'sp3',
        'CG'  => 'sp2',
        'CD1' => 'sp2',
        'CD2' => 'sp2',
        'NE1' => 'sp2',
        'CE2' => 'sp2',
        'CE3' => 'sp2',
        'CZ2' => 'sp2',
        'CZ3' => 'sp2',
        'CH2' => 'sp2',
        'OXT' => 'sp2',
    },
    'THR' => {
        'CA'  => 'sp3',
        'C'   => 'sp2',
        'O'   => 'sp2',
        'CB'  => 'sp3',
        'OG1' => 'sp3',
        'CG2' => 'sp3',
        'OXT' => 'sp2',
    },
    'ASN' => {
        'CA'  => 'sp3',
        'C'   => 'sp2',
        'O'   => 'sp2',
        'CB'  => 'sp3',
        'CG'  => 'sp2',
        'OD1' => 'sp2',
        'ND2' => 'sp2',
        'OXT' => 'sp2',
    },
    'GLN' => {
        'CA'  => 'sp3',
        'C'   => 'sp2',
        'O'   => 'sp2',
        'CB'  => 'sp3',
        'CG'  => 'sp3',
        'CD'  => 'sp2',
        'OE1' => 'sp2',
        'NE2' => 'sp2',
        'OXT' => 'sp2',
    },
);

our @ROTATABLE_RESIDUE_NAMES =
    qw(XAA SER ARG HIS LYS ASP GLU CYS VAL ILE LEU MET PHE TYR TRP THR ASN GLN);

# ------------------------ Molecule-related functions ------------------------- #

#
# Creates a hash where for each unique residue key proper atom ids are assigned.
# Input:
#     $atom_site - $atom_site - atom data structure.
# Output:
#     %residue_atom_ids - hash of unique residue key and corresponding atom ids.
#

sub identify_residue_atoms
{
    my ( $atom_site, $options ) = @_;
    my ( $check_atom_names ) = $options->{'check_atom_names'};

    my $split_groups = split_by( { 'atom_site' => $atom_site,
                                   'attributes' => [ 'label_seq_id',
                                                     'label_asym_id',
                                                     'pdbx_PDB_model_num', ] } );

    my %residue_atom_ids;
    for my $atom_id ( keys %{ $atom_site } ) {
        my $unique_residue_key = unique_residue_key( $atom_site->{$atom_id} );
        my ( $residue_id, $residue_chain, $pdbx_model, $alt_id ) =
            split /,/, $unique_residue_key;
        my $split_group_entry = "${residue_id},${residue_chain},${pdbx_model}";
        my $related_residue_atom_ids = $split_groups->{$split_group_entry};

        # Splits related residues into alt id groups. Atom ids can be redundant.
        for my $related_atom_id ( @{ $related_residue_atom_ids } ) {
            my $related_alt_id =
                $atom_site->{$related_atom_id}{'label_alt_id'};

            if( $atom_id ne $related_atom_id ) {
                next if( $check_atom_names &&
                         $atom_site->{$atom_id}{'label_atom_id'} eq
                         $atom_site->{$related_atom_id}{'label_atom_id'} );

                if( $alt_id eq '.' ) {
                    push @{ $residue_atom_ids{$unique_residue_key} },
                        $related_atom_id;
                } elsif( $alt_id eq $related_alt_id || $related_alt_id eq '.') {
                    push @{ $residue_atom_ids{$unique_residue_key} },
                        $related_atom_id;
                }
            }
        }
    }

    return \%residue_atom_ids;
}

1;
