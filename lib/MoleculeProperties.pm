package MoleculeProperties;

use strict;
use warnings;

use Exporter qw( import );
our @EXPORT_OK = qw( %HYBRIDIZATION
                     %HYDROGEN_NAMES
                     %ROTATABLE_BONDS );

# ------------------------------- Molecule properties ------------------------- #

#
# Stores molecule properties, such as rotatable bonds of certain sidechains.
#

our %ROTATABLE_BONDS = (
    'SER' => {
	'HB2'  => [ [ 'CA', 'CB' ], [ 'CB', 'HB2' ] ],
	'HB3'  => [ [ 'CA', 'CB' ], [ 'CB', 'HB3' ] ],
	'OG'   => [ [ 'CA', 'CB' ], [ 'CB', 'OG' ] ],
	'HG'   => [ [ 'CA', 'CB' ], [ 'CB', 'OG' ], [ 'OG', 'HG' ] ]
    },
    'PHE' => {
	'CD2' => [ [ 'CA', 'CB' ], [ 'CB', 'CG' ], [ 'CG', 'CD1' ] ],
	'CE2' => [ [ 'CA', 'CB' ], [ 'CB', 'CG' ], [ 'CG', 'CD1' ] ],
	'CD1' => [ [ 'CA', 'CB' ], [ 'CB', 'CG' ], [ 'CG', 'CD1' ] ],
	'CZ'  => [ [ 'CA', 'CB' ], [ 'CB', 'CG' ], [ 'CG', 'CD1' ] ],
	'CE1' => [ [ 'CA', 'CB' ], [ 'CB', 'CG' ], [ 'CG', 'CD1' ] ],
	'CG'  => [ [ 'CA', 'CB' ], [ 'CB', 'CG' ] ]
    },
    'TRP' => {
	'CE3' => [ [ 'CA', 'CB' ], [ 'CB', 'CG' ], [ 'CG', 'CD1' ] ],
	'CH2' => [ [ 'CA', 'CB' ], [ 'CB', 'CG' ], [ 'CG', 'CD1' ] ],
	'CE2' => [ [ 'CA', 'CB' ], [ 'CB', 'CG' ], [ 'CG', 'CD1' ] ],
	'CD2' => [ [ 'CA', 'CB' ], [ 'CB', 'CG' ], [ 'CG', 'CD1' ] ],
	'CD1' => [ [ 'CA', 'CB' ], [ 'CB', 'CG' ], [ 'CG', 'CD1' ] ],
	'CZ2' => [ [ 'CA', 'CB' ], [ 'CB', 'CG' ], [ 'CG', 'CD1' ] ],
	'CZ3' => [ [ 'CA', 'CB' ], [ 'CB', 'CG' ], [ 'CG', 'CD1' ] ],
	'NE1' => [ [ 'CA', 'CB' ], [ 'CB', 'CG' ], [ 'CG', 'CD1' ] ],
	'CG'  => [ [ 'CA', 'CB' ], [ 'CB', 'CG' ] ]
    },
    'TYR' => {
	'CD1' => [ [ 'CA', 'CB' ], [ 'CB', 'CG' ], [ 'CG', 'CD1' ] ],
	'CD2' => [ [ 'CA', 'CB' ], [ 'CB', 'CG' ], [ 'CG', 'CD1' ] ],
	'CE2' => [ [ 'CA', 'CB' ], [ 'CB', 'CG' ], [ 'CG', 'CD1' ] ],
	'CZ'  => [ [ 'CA', 'CB' ], [ 'CB', 'CG' ], [ 'CD1' ] ],
	'CG'  => [ [ 'CA', 'CB' ], [ 'CB', 'CG' ] ],
	'CE1' => [ [ 'CA', 'CB' ], [ 'CB', 'CG' ], [ 'CG', 'CD1' ] ],
	'OH'  => [ [ 'CA', 'CB' ], [ 'CB', 'CG' ], [ 'CG', 'CD1' ] ]
    },
    'THR' => {
	'CG2' => [ [ 'CA', 'CB' ], [ 'CB', 'OG1' ] ],
	'OG1' => [ [ 'CA', 'CB' ], [ 'CB', 'OG1' ] ]
    },
    'ILE' => {
	'CD1' => [ [ 'CA', 'CB' ], [ 'CB', 'CG1' ], [ 'CG1', 'CD1' ] ],
	'CG1' => [ [ 'CA', 'CB' ], [ 'CB', 'CG1' ] ],
	'CG2' => [ [ 'CA', 'CB' ], [ 'CB', 'CG1' ] ]
    },
    'ARG' => {
	'CG'  => [ [ 'CA', 'CB' ], [ 'CB', 'CG' ] ],
	'CZ'  => [ [ 'CA', 'CB' ], [ 'CB', 'CG' ], [ 'CG', 'CD' ], [ 'CD', 'NE' ], [ 'NE', 'CZ' ] ],
	'NE'  => [ [ 'CA', 'CB' ], [ 'CB', 'CG' ], [ 'CG', 'CD' ], [ 'CD', 'NE' ] ],
	'CD'  => [ [ 'CA', 'CB' ], [ 'CB', 'CG' ], [ 'CG', 'CD' ] ],
	'NH1' => [ [ 'CA', 'CB' ], [ 'CB', 'CG' ], [ 'CG', 'CD' ], [ 'CD', 'NE' ], [ 'NE', 'CZ' ] ],
	'NH2' => [ [ 'CA', 'CB' ], [ 'CB', 'CG' ], [ 'CG', 'CD' ], [ 'CD', 'NE' ], [ 'NE', 'CZ' ] ]
    },
    'MET' => {
	'CE'  => [ [ 'CA', 'CB' ], [ 'CB', 'CG' ], [ 'CG', 'SD' ], [ 'SD', 'CE' ] ],
	'SD'  => [ [ 'CA', 'CB' ], [ 'CB', 'CG' ], [ 'CG', 'SD' ] ],
	'CG'  => [ [ 'CA', 'CB' ], [ 'CB', 'CG' ] ]
    },
    'ASN' => {
	'CG'  => [ [ 'CA', 'CB' ], [ 'CB', 'CG' ] ],
	'ND2' => [ [ 'CA', 'CB' ], [ 'CB', 'CG' ], [ 'CG', 'OD1' ] ],
	'OD1' => [ [ 'CA', 'CB' ], [ 'CB', 'CG' ], [ 'CG', 'OD1' ] ]
    },
    'LEU' => {
	'CD2' => [ [ 'CA', 'CB' ], [ 'CB', 'CG' ], [ 'CG', 'CD1' ] ],
	'CD1' => [ [ 'CA', 'CB' ], [ 'CB', 'CG' ], [ 'CG', 'CD1' ] ],
	'CG'  => [ [ 'CA', 'CB' ], [ 'CB', 'CG' ] ]
    },
    'GLN' => {
	'NE2' => [ [ 'CA', 'CB' ], [ 'CB', 'CG' ], [ 'CG', 'CD' ], [ 'CD', 'OE1' ] ],
	'OE1' => [ [ 'CA', 'CB' ], [ 'CB', 'CG' ], [ 'CG', 'CD' ], [ 'CD', 'OE1' ] ],
	'CD'  => [ [ 'CA', 'CB' ], [ 'CB', 'CG' ], [ 'CG', 'CD' ] ],
	'CG'  => [ [ 'CA', 'CB' ], [ 'CB', 'CG' ] ]
    },
    'CYS' => {
	'SG'  => [ [ 'CA', 'CB' ], [ 'CB', 'SG' ] ]
    },
    'ASP' => {
	'OD1' => [ [ 'CA', 'CB' ], [ 'CB', 'CG' ], [ 'CG', 'OD1' ] ],
	'CG'  => [ [ 'CA', 'CB' ], [ 'CB', 'CG' ] ],
	'OD2' => [ [ 'CA', 'CB' ], [ 'CB', 'CG' ], [ 'CG', 'OD1' ] ]
    },
    'GLU' => {
	'CG'  => [ [ 'CA', 'CB' ], [ 'CB', 'CG' ] ],
	'OE2' => [ [ 'CA', 'CB' ], [ 'CB', 'CG' ], [ 'CG', 'CD' ], [ 'CD', 'OE1' ] ],
	'CD'  => [ [ 'CA', 'CB' ], [ 'CB', 'CG' ], [ 'CG', 'CD' ] ],
	'OE1' => [ [ 'CA', 'CB' ], [ 'CB', 'CG' ], [ 'CG', 'CD' ], [ 'CD', 'OE1' ] ]
    },
    'VAL' => {
	'CG2' => [ [ 'CA', 'CB' ], [ 'CB', 'CG1' ] ],
	'CG1' => [ [ 'CA', 'CB' ], [ 'CB', 'CG1' ] ]
    },
    'HIS' => {
	'ND1' => [ [ 'CA', 'CB' ], [ 'CB', 'CG' ], [ 'CG', 'ND1' ] ],
	'CE1' => [ [ 'CA', 'CB' ], [ 'CB', 'CG' ], [ 'CG', 'ND1' ] ],
	'CG'  => [ [ 'CA', 'CB' ], [ 'CB', 'CG' ] ],
	'NE2' => [ [ 'CA', 'CB' ], [ 'CB', 'CG' ], [ 'CG', 'ND1' ] ],
	'CD2' => [ [ 'CA', 'CB' ], [ 'CB', 'CG' ], [ 'CG', 'ND1' ] ]
    },
    'LYS' => {
	'CG'  => [ [ 'CA', 'CB' ], [ 'CB', 'CG' ] ],
	'CD'  => [ [ 'CA', 'CB' ], [ 'CB', 'CG' ], [ 'CG', 'CD' ] ],
	'CE'  => [ [ 'CA', 'CB' ], [ 'CB', 'CG' ], [ 'CG', 'CD' ], [ 'CD', 'CE' ] ],
	'NZ'  => [ [ 'CA', 'CB' ], [ 'CB', 'CG' ], [ 'CG', 'CD' ], [ 'CD', 'CE' ], [ 'CE', 'NZ' ] ]
    } );

our %HYBRIDIZATION = (
    'SER' => {
	'N'  => 'sp3',
	'CA' => 'sp3',
	'C'  => 'sp2',
	'O'  => 'sp2',
	'OXT'=> 'sp3',
	'CB' => 'sp3',
	'OG' => 'sp3'
    },
    'ARG' => {
	'N'   => 'sp3',
	'CA'  => 'sp3',
	'C'   => 'sp2',
	'O'   => 'sp2',
	'OXT' => 'sp3',
	'CB'  => 'sp3',
	'CG'  => 'sp3',
	'CD'  => 'sp3',
	'NE'  => 'sp2',
	'CZ'  => 'sp2',
	'NH1' => 'sp2',
	'NH2' => 'sp2'
    },
    'HIS' => {
	'N'   => 'sp3',
	'CA'  => 'sp3',
	'C'   => 'sp2',
	'O'   => 'sp2',
	'OXT' => 'sp3',
	'CB'  => 'sp3',
	'CG'  => 'sp2',
	'ND1' => 'sp2',
	'CD2' => 'sp2',
	'CE1' => 'sp2',
	'NE2' => 'sp2'
    },
    'LYS' => {
	'N'   => 'sp3',
	'CA'  => 'sp3',
	'C'   => 'sp2',
	'O'   => 'sp2',
	'OXT' => 'sp3',
	'CB'  => 'sp3',
	'CD'  => 'sp3',
	'CG'  => 'sp3',
	'CE'  => 'sp3',
	'NZ'  => 'sp3'
    },
    'ASP' => {
	'N'   => 'sp3',
	'CA'  => 'sp3',
	'C'   => 'sp2',
	'O'   => 'sp2',
	'OXT' => 'sp3',
	'CB'  => 'sp3',
	'CG'  => 'sp2',
	'OD1' => 'sp2',
	'OD2' => 'sp2'
    },
    'GLU' => {
	'N'   => 'sp3',
	'CA'  => 'sp3',
	'C'   => 'sp2',
	'O'   => 'sp2',
	'OXT' => 'sp3',
	'CB'  => 'sp3',
	'CG'  => 'sp3',
	'CD'  => 'sp2',
	'OE1' => 'sp2',
	'OE2' => 'sp2'
    },
    'CYS' => {
	'N'   => 'sp3',
	'CA'  => 'sp3',
	'C'   => 'sp2',
	'O'   => 'sp2',
	'OXT' => 'sp3',
	'CB'  => 'sp3',
	'SG'  => 'sp3'
    },
    'GLY' => {
	'N'   => 'sp3',
	'CA'  => 'sp3',
	'C'   => 'sp2',
	'O'   => 'sp2',
	'OXT' => 'sp3'
    },
    'PRO' => {
	'N'   => 'sp3',
	'CA'  => 'sp3',
	'C'   => 'sp2',
	'O'   => 'sp2',
	'OXT' => 'sp3',
	'CB'  => 'sp3',
	'CG'  => 'sp3',
	'CD'  => 'sp3'
    },
    'ALA' => {
	'N'   => 'sp3',
	'CA'  => 'sp3',
	'C'   => 'sp2',
	'O'   => 'sp2',
	'OXT' => 'sp3',
	'CB'  => 'sp3',
    },
    'VAL' => {
	'N'   => 'sp3',
	'CA'  => 'sp3',
	'C'   => 'sp2',
	'O'   => 'sp2',
	'OXT' => 'sp3',
	'CB'  => 'sp3',
	'CG1' => 'sp3',
	'CG2' => 'sp3'
    },
    'ILE' => {
	'N'   => 'sp3',
	'CA'  => 'sp3',
	'C'   => 'sp2',
	'O'   => 'sp2',
	'OXT '=> 'sp3',
	'CB'  => 'sp3',
	'CG1' => 'sp3',
	'CG2' => 'sp3',
	'CD1' => 'sp3'
    },
    'LEU' => {
	'N'   => 'sp3',
	'CA'  => 'sp3',
	'C'   => 'sp2',
	'O'   => 'sp2',
	'OXT '=> 'sp3',
	'CB'  => 'sp3',
	'CG'  => 'sp3',
	'CD1' => 'sp3',
	'CD2' => 'sp3'
    },
    'MET' => {
	'N'   => 'sp3',
	'CA'  => 'sp3',
	'C'   => 'sp2',
	'O'   => 'sp2',
	'OXT' => 'sp3',
	'CB'  => 'sp3',
	'CG'  => 'sp3',
	'SD'  => 'sp3',
	'CE'  => 'sp3'
    },
    'PHE' => {
	'N'   => 'sp3',
	'CA'  => 'sp3',
	'C'   => 'sp2',
	'O'   => 'sp2',
	'OXT' => 'sp3',
	'CB'  => 'sp3',
	'CG'  => 'sp2',
	'CD1' => 'sp2',
	'CD2' => 'sp2',
	'CE1' => 'sp2',
	'CE2' => 'sp2',
	'CZ'  => 'sp2',
    },
    'TYR' => {
	'N'   => 'sp3',
	'CA'  => 'sp3',
	'C'   => 'sp2',
	'O'   => 'sp2',
	'OXT' => 'sp3',
	'CB'  => 'sp3',
	'CG'  => 'sp2',
	'CE1' => 'sp2',
	'CE2' => 'sp2',
	'CD1' => 'sp2',
	'CD2' => 'sp2',
	'CZ'  => 'sp2',
	'OH'  => 'sp3'
    },
    'TRP' => {
	'N'   => 'sp3',
	'CA'  => 'sp3',
	'C'   => 'sp2',
	'O'   => 'sp2',
	'OXT' => 'sp3',
	'CB'  => 'sp3',
	'CG'  => 'sp2',
	'CD1' => 'sp2',
	'CD2' => 'sp2',
	'NE1' => 'sp2',
	'CE2' => 'sp2',
	'CE3' => 'sp2',
	'CZ2' => 'sp2',
	'CZ3' => 'sp2',
	'CH2' => 'sp2'
    },
    'THR' => {
	'N'   => 'sp3',
	'CA'  => 'sp3',
	'C'   => 'sp2',
	'O'   => 'sp2',
	'OXT' => 'sp3',
	'CB'  => 'sp3',
	'OG1' => 'sp3',
	'CG2' => 'sp3'
    },
    'ASN' => {
	'N'   => 'sp3',
	'CA'  => 'sp3',
	'C'   => 'sp2',
	'O'   => 'sp2',
	'OXT' => 'sp3',
	'CB'  => 'sp3',
	'CG'  => 'sp2',
	'OD1' => 'sp2',
	'ND2' => 'sp2'
    },
    'GLN' => {
	'N'   => 'sp3',
	'CA'  => 'sp3',
	'C'   => 'sp2',
	'O'   => 'sp2',
	'OXT' => 'sp3',
	'CB'  => 'sp3',
	'CG'  => 'sp3',
	'CD'  => 'sp2',
	'OE1' => 'sp2',
	'NE2' => 'sp2'
    },
);

our %HYDROGEN_NAMES = (
    'SER' => {
	'N'   => [ 'H', 'H2' ],
	'CA'  => [ 'HA' ],
	'OXT' => [ 'HXT' ],
	'CB'  => [ 'HB2', 'HB3' ],
	'OG'  => [ 'HG' ]
    },
    'ARG' => {
	'N'   => [ 'H', 'H2' ],
	'CA'  => [ 'HA' ],
	'OXT' => [ 'HXT' ],
	'CB'  => [ 'HB2', 'HB3' ],
	'CG'  => [ 'HG2', 'HG3' ],
	'CD'  => [ 'HD2', 'HD3' ],
	'NE'  => [ 'HE' ],
	'NH1' => [ 'HH11', 'HH12' ],
	'NH2' => [ 'HH21', 'HH22' ]
    },
    'HIS' => {
	'N'   => [ 'H', 'H2' ],
	'CA'  => [ 'HA' ],
	'OXT' => [ 'HXT' ],
	'CB'  => [ 'HB2', 'HB3' ],
	'ND1' => [ 'HD1' ],
	'CD2' => [ 'HD2' ],
	'CE1' => [ 'HE1' ],
	'NE2' => [ 'HE2' ]
    },
    'LYS' => {
	'N'   => [ 'H', 'H2' ],
	'CA'  => [ 'HA' ],
	'OXT' => [ 'HXT' ],
	'CB'  => [ 'HB2', 'HB3' ],
	'CG'  => [ 'HG2', 'HG3' ],
	'CD'  => [ 'HD2', 'HD3' ],
	'CE'  => [ 'HE2', 'HE3' ],
	'NZ'  => [ 'HZ1', 'HZ2' ]
    },
    'ASP' => {
	'N'   => [ 'H', 'H2' ],
	'CA'  => [ 'HA' ],
	'OXT' => [ 'HXT' ],
	'CB'  => [ 'HB2', 'HB3' ],
	'OD2' => [ 'HD2' ]
    },
    'GLU' => {
	'N'   => [ 'H', 'H2' ],
	'CA'  => [ 'HA' ],
	'OXT' => [ 'HXT' ],
	'CB'  => [ 'HB2', 'HB3' ],
	'CG'  => [ 'HG2', 'HG3' ],
	'OE2' => [ 'HE2' ]
    },
    'CYS' => {
	'N'   => [ 'H', 'H2' ],
	'CA'  => [ 'HA' ],
	'OXT' => [ 'HXT' ],
	'CB'  => [ 'HB2', 'HB3' ],
	'SG'  => [ 'HG' ]
    },
    'GLY' => {
	'N'   => [ 'H', 'H2' ],
	'CA'  => [ 'HA2', 'HA3' ],
	'OXT' => [ 'HXT' ]
    },
    'PRO' => {
	'N'   => [ 'H' ],
	'CA'  => [ 'HA' ],
	'CB'  => [ 'HB2', 'HB3' ],
	'CG'  => [ 'HG2', 'HG3' ],
	'CD'  => [ 'HD2', 'HD3' ]
    },
    'ALA' => {
	'N'   => [ 'H', 'H2' ],
	'CA'  => [ 'HA' ],
	'OXT' => [ 'HXT' ],
	'CB'  => [ 'HB1', 'HB2', 'HB3' ]
    },
    'VAL' => {
	'N'   => [ 'H', 'H2' ],
	'CA'  => [ 'HA' ],
	'OXT' => [ 'HXT' ],
	'CB'  => [ 'HB' ],
	'CG1' => [ 'HG11', 'HG12', 'HG13' ],
	'CG2' => [ 'HG21', 'HG22', 'HG23' ]
    },
    'ILE' => {
	'N'   => [ 'H', 'H2' ],
	'CA'  => [ 'HA' ],
	'OXT' => [ 'HXT' ],
	'CB'  => [ 'HB' ],
	'CG1' => [ 'HG12', 'HG13' ],
	'CG2' => [ 'HG21', 'HG22', 'HG23' ],
	'CD1' => [ 'HD11', 'HD12', 'HD13' ]
    },
    'LEU' => {
	'N'   => [ 'H', 'H2' ],
	'CA'  => [ 'HA' ],
	'OXT' => [ 'HXT' ],
	'CB'  => [ 'HB2', 'HB3' ],
	'CG'  => [ 'HG' ],
	'CD1' => [ 'HD11', 'HD12', 'HD13' ],
	'CD2' => [ 'HD21', 'HD22', 'HD23' ]
    },
    'MET' => {
	'N'   => [ 'H', 'H2' ],
	'CA'  => [ 'HA' ],
	'OXT' => [ 'HXT' ],
	'CB'  => [ 'HB2', 'HB3' ],
	'CG'  => [ 'HG2', 'HG3' ],
	'CE'  => [ 'HE1', 'HE2', 'HE3' ]
    },
    'PHE' => {
	'N'   => [ 'H', 'H2' ],
	'CA'  => [ 'HA' ],
	'OXT' => [ 'HXT' ],
	'CB'  => [ 'HB2', 'HB3' ],
	'CD1' => [ 'HD1' ],
	'CD2' => [ 'HD2' ],
	'CE1' => [ 'HE1' ],
	'CE2' => [ 'HE2' ],
	'CZ'  => [ 'HZ' ]
    },
    'TYR' => {
	'N'   => [ 'H', 'H2' ],
	'CA'  => [ 'HA' ],
	'OXT' => [ 'HXT' ],
	'CB'  => [ 'HB2', 'HB3' ],
	'CD1' => [ 'HD1' ],
	'CD2' => [ 'HD2' ],
	'CE1' => [ 'HE1' ],
	'CE2' => [ 'HE2' ],
	'OH'  => [ 'HH' ]
    },
    'TRP' => {
	'N'   => [ 'H', 'H2' ],
	'CA'  => [ 'HA' ],
	'OXT' => [ 'HXT' ],
	'CB'  => [ 'HB2', 'HB3' ],
	'CD1' => [ 'HD1' ],
	'NE1' => [ 'HE1' ],
	'CE3' => [ 'HE3' ],
	'CZ2' => [ 'HZ2' ],
	'CZ3' => [ 'HZ3' ],
	'CH2' => [ 'HH2' ]
    },
    'THR' => {
	'N'   => [ 'H', 'H2' ],
	'CA'  => [ 'HA' ],
	'OXT' => [ 'HXT' ],
	'CB'  => [ 'HB' ],
	'OG1' => [ 'HG1' ],
	'CG2' => [ 'HG21', 'HG22', 'HG23' ]
    },
    'ASN' => {
	'N'   => [ 'H', 'H2' ],
	'CA'  => [ 'HA' ],
	'OXT' => [ 'HXT' ],
	'CB'  => [ 'HB2', 'HB3' ],
	'ND2' => [ 'HD21', 'HD22' ]
    },
    'GLN' => {
	'N'   => [ 'H', 'H2' ],
	'CA'  => [ 'HA' ],
	'OXT' => [ 'HXT' ],
	'CB'  => [ 'HB2', 'HB3' ],
	'CG'  => [ 'HG2', 'HG3' ],
	'NE2' => [ 'HE21', 'HE22' ]
    }
);

1;
