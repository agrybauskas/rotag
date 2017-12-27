package MoleculeProperties;

use strict;
use warnings;

use Exporter qw( import );
our @EXPORT_OK = qw( %BOND_TYPES
                     %HYBRIDIZATION
                     %HYDROGEN_NAMES
                     %ROTATABLE_BONDS );

use AtomProperties qw( %ATOMS );

# ------------------------------- Molecule properties ------------------------- #

#
# Stores molecule properties, such as rotatable bonds of certain sidechains,
# hybridization angles, hydrogen names for specific residues.
#

our %ROTATABLE_BONDS = (
    'SER' => {
	'HB2' => [ [ 'CA', 'CB' ], [ 'CB', 'OG' ] ],
	'HB3' => [ [ 'CA', 'CB' ], [ 'CB', 'OG' ] ],
	'OG'  => [ [ 'CA', 'CB' ], [ 'CB', 'OG' ] ],
	'HG'  => [ [ 'CA', 'CB' ], [ 'CB', 'OG' ], [ 'OG', 'HG' ] ]
    },
    'ARG' => {
	'HB2' => [ [ 'CA', 'CB' ], [ 'CB', 'OG' ] ],
	'HB3' => [ [ 'CA', 'CB' ], [ 'CB', 'OG' ] ],
	'CG'  => [ [ 'CA', 'CB' ], [ 'CB', 'CG' ] ],
	'HG2' => [ [ 'CA', 'CB' ], [ 'CB', 'CG' ], [ 'CG', 'CD' ] ],
	'HG3' => [ [ 'CA', 'CB' ], [ 'CB', 'CG' ], [ 'CG', 'CD' ] ],
	'CD'  => [ [ 'CA', 'CB' ], [ 'CB', 'CG' ], [ 'CG', 'CD' ] ],
	'HD2' => [ [ 'CA', 'CB' ], [ 'CB', 'CG' ], [ 'CG', 'CD' ], [ 'CD', 'NE' ] ],
	'HD3' => [ [ 'CA', 'CB' ], [ 'CB', 'CG' ], [ 'CG', 'CD' ], [ 'CD', 'NE' ] ],
	'NE'  => [ [ 'CA', 'CB' ], [ 'CB', 'CG' ], [ 'CG', 'CD' ], [ 'CD', 'NE' ] ],
	'HE'  => [ [ 'CA', 'CB' ], [ 'CB', 'CG' ], [ 'CG', 'CD' ], [ 'CD', 'NE' ], [ 'NE', 'CZ' ] ],
	'CZ'  => [ [ 'CA', 'CB' ], [ 'CB', 'CG' ], [ 'CG', 'CD' ], [ 'CD', 'NE' ], [ 'NE', 'CZ' ] ],
	'NH1' => [ [ 'CA', 'CB' ], [ 'CB', 'CG' ], [ 'CG', 'CD' ], [ 'CD', 'NE' ], [ 'NE', 'CZ' ] ],
	'HH11'=> [ [ 'CA', 'CB' ], [ 'CB', 'CG' ], [ 'CG', 'CD' ], [ 'CD', 'NE' ], [ 'NE', 'CZ' ] ],
	'HH12'=> [ [ 'CA', 'CB' ], [ 'CB', 'CG' ], [ 'CG', 'CD' ], [ 'CD', 'NE' ], [ 'NE', 'CZ' ] ],
	'NH2' => [ [ 'CA', 'CB' ], [ 'CB', 'CG' ], [ 'CG', 'CD' ], [ 'CD', 'NE' ], [ 'NE', 'CZ' ] ],
	'HH21'=> [ [ 'CA', 'CB' ], [ 'CB', 'CG' ], [ 'CG', 'CD' ], [ 'CD', 'NE' ], [ 'NE', 'CZ' ] ]
    },
    'HIS' => {
	'HB2' => [ [ 'CA', 'CB' ], [ 'CB', 'CG' ] ],
	'HB3' => [ [ 'CA', 'CB' ], [ 'CB', 'CG' ] ],
	'CG'  => [ [ 'CA', 'CB' ], [ 'CB', 'CG' ] ],
	'HD1' => [ [ 'CA', 'CB' ], [ 'CB', 'CG' ], [ 'CG', 'ND1' ] ],
	'ND1' => [ [ 'CA', 'CB' ], [ 'CB', 'CG' ], [ 'CG', 'ND1' ] ],
	'HD2' => [ [ 'CA', 'CB' ], [ 'CB', 'CG' ], [ 'CG', 'ND1' ] ],
	'CD2' => [ [ 'CA', 'CB' ], [ 'CB', 'CG' ], [ 'CG', 'ND1' ] ],
	'HE1' => [ [ 'CA', 'CB' ], [ 'CB', 'CG' ], [ 'CG', 'ND1' ] ],
	'CE1' => [ [ 'CA', 'CB' ], [ 'CB', 'CG' ], [ 'CG', 'ND1' ] ],
	'HE2' => [ [ 'CA', 'CB' ], [ 'CB', 'CG' ], [ 'CG', 'ND1' ] ],
	'NE2' => [ [ 'CA', 'CB' ], [ 'CB', 'CG' ], [ 'CG', 'ND1' ] ]
    },
    'LYS' => {
	'HB2' => [ [ 'CA', 'CB' ], [ 'CB', 'CG' ] ],
	'HB3' => [ [ 'CA', 'CB' ], [ 'CB', 'CG' ] ],
	'CG'  => [ [ 'CA', 'CB' ], [ 'CB', 'CG' ] ],
	'HG2' => [ [ 'CA', 'CB' ], [ 'CB', 'CG' ], [ 'CG', 'CD' ] ],
	'HG3' => [ [ 'CA', 'CB' ], [ 'CB', 'CG' ], [ 'CG', 'CD' ] ],
	'CD'  => [ [ 'CA', 'CB' ], [ 'CB', 'CG' ], [ 'CG', 'CD' ] ],
	'HD2' => [ [ 'CA', 'CB' ], [ 'CB', 'CG' ], [ 'CG', 'CD' ], [ 'CD', 'CE' ] ],
	'HD3' => [ [ 'CA', 'CB' ], [ 'CB', 'CG' ], [ 'CG', 'CD' ], [ 'CD', 'CE' ] ],
	'CE'  => [ [ 'CA', 'CB' ], [ 'CB', 'CG' ], [ 'CG', 'CD' ], [ 'CD', 'CE' ] ],
	'HE2' => [ [ 'CA', 'CB' ], [ 'CB', 'CG' ], [ 'CG', 'CD' ], [ 'CD', 'CE' ], [ 'CE', 'NZ' ] ],
	'HE3' => [ [ 'CA', 'CB' ], [ 'CB', 'CG' ], [ 'CG', 'CD' ], [ 'CD', 'CE' ], [ 'CE', 'NZ' ] ],
	'NZ'  => [ [ 'CA', 'CB' ], [ 'CB', 'CG' ], [ 'CG', 'CD' ], [ 'CD', 'CE' ], [ 'CE', 'NZ' ] ],
	'HZ1' => [ [ 'CA', 'CB' ], [ 'CB', 'CG' ], [ 'CG', 'CD' ], [ 'CD', 'CE' ], [ 'CE', 'NZ' ], [ 'NZ', 'HZ1' ] ],
	'HZ2' => [ [ 'CA', 'CB' ], [ 'CB', 'CG' ], [ 'CG', 'CD' ], [ 'CD', 'CE' ], [ 'CE', 'NZ' ], [ 'NZ', 'HZ1' ] ]
    },
    'ASP' => {
	'HB2' => [ [ 'CA', 'CB' ], [ 'CB', 'CG' ] ],
	'HB3' => [ [ 'CA', 'CB' ], [ 'CB', 'CG' ] ],
	'CG'  => [ [ 'CA', 'CB' ], [ 'CB', 'CG' ] ],
	'OD1' => [ [ 'CA', 'CB' ], [ 'CB', 'CG' ], [ 'CG', 'OD1' ] ],
	'OD2' => [ [ 'CA', 'CB' ], [ 'CB', 'CG' ], [ 'CG', 'OD1' ] ],
	'HD2' => [ [ 'CA', 'CB' ], [ 'CB', 'CG' ], [ 'CG', 'OD1' ] ]
    },
    'GLU' => {
	'HB2' => [ [ 'CA', 'CB' ], [ 'CB', 'CG' ] ],
	'HB3' => [ [ 'CA', 'CB' ], [ 'CB', 'CG' ] ],
	'CG'  => [ [ 'CA', 'CB' ], [ 'CB', 'CG' ] ],
	'HG2' => [ [ 'CA', 'CB' ], [ 'CB', 'CG' ], [ 'CG', 'CD' ] ],
	'HG3' => [ [ 'CA', 'CB' ], [ 'CB', 'CG' ], [ 'CG', 'CD' ] ],
	'CD'  => [ [ 'CA', 'CB' ], [ 'CB', 'CG' ], [ 'CG', 'CD' ] ],
	'OE1' => [ [ 'CA', 'CB' ], [ 'CB', 'CG' ], [ 'CG', 'CD' ], [ 'CD', 'OE1' ] ],
	'OE2' => [ [ 'CA', 'CB' ], [ 'CB', 'CG' ], [ 'CG', 'CD' ], [ 'CD', 'OE1' ] ],
	'HE2' => [ [ 'CA', 'CB' ], [ 'CB', 'CG' ], [ 'CG', 'CD' ], [ 'CD', 'OE1' ] ]
    },
    'CYS' => {
	'HB2' => [ [ 'CA', 'CB' ], [ 'CB', 'SG' ] ],
	'HB3' => [ [ 'CA', 'CB' ], [ 'CB', 'SG' ] ],
	'SG'  => [ [ 'CA', 'CB' ], [ 'CB', 'SG' ] ],
	'HG'  => [ [ 'CA', 'CB' ], [ 'CB', 'SG' ], [ 'SG', 'HG' ] ]
    },
    'VAL' => {
	'HB'  => [ [ 'CA', 'CB' ], [ 'CB', 'CG1' ] ],
	'CG1' => [ [ 'CA', 'CB' ], [ 'CB', 'CG1' ] ],
	'HG11'=> [ [ 'CA', 'CB' ], [ 'CB', 'CG1' ], [ 'CG1', 'HG11' ] ],
	'HG12'=> [ [ 'CA', 'CB' ], [ 'CB', 'CG1' ], [ 'CG1', 'HG11' ] ],
	'HG13'=> [ [ 'CA', 'CB' ], [ 'CB', 'CG1' ], [ 'CG1', 'HG11' ] ],
	'CG2' => [ [ 'CA', 'CB' ], [ 'CB', 'CG1' ] ],
	'HG21'=> [ [ 'CA', 'CB' ], [ 'CB', 'CG2' ], [ 'CG2', 'HG21' ] ],
	'HG22'=> [ [ 'CA', 'CB' ], [ 'CB', 'CG2' ], [ 'CG2', 'HG21' ] ],
	'HG23'=> [ [ 'CA', 'CB' ], [ 'CB', 'CG2' ], [ 'CG2', 'HG21' ] ]
    },
    'ILE' => {
	'HB'  => [ [ 'CA', 'CB' ], [ 'CB', 'CG1' ] ],
	'CG1' => [ [ 'CA', 'CB' ], [ 'CB', 'CG1' ] ],
	'HG11'=> [ [ 'CA', 'CB' ], [ 'CB', 'CG1' ], [ 'CG1', 'CD1' ] ],
	'HG12'=> [ [ 'CA', 'CB' ], [ 'CB', 'CG1' ], [ 'CG1', 'CD1' ] ],
	'HG13'=> [ [ 'CA', 'CB' ], [ 'CB', 'CG1' ], [ 'CG1', 'CD1' ] ],
	'CD1' => [ [ 'CA', 'CB' ], [ 'CB', 'CG1' ], [ 'CG1', 'CD1' ] ],
	'HD11'=> [ [ 'CA', 'CB' ], [ 'CB', 'CG1' ], [ 'CG1', 'CD1' ], [ 'CD1', 'HD11' ] ],
	'HD12'=> [ [ 'CA', 'CB' ], [ 'CB', 'CG1' ], [ 'CG1', 'CD1' ], [ 'CD1', 'HD11' ] ],
	'HD13'=> [ [ 'CA', 'CB' ], [ 'CB', 'CG1' ], [ 'CG1', 'CD1' ], [ 'CD1', 'HD11' ] ],
	'CG2' => [ [ 'CA', 'CB' ], [ 'CB', 'CG1' ] ],
	'HG21'=> [ [ 'CA', 'CB' ], [ 'CB', 'CG2' ], [ 'CG2', 'HG21' ] ],
	'HG22'=> [ [ 'CA', 'CB' ], [ 'CB', 'CG2' ], [ 'CG2', 'HG21' ] ],
	'HG23'=> [ [ 'CA', 'CB' ], [ 'CB', 'CG2' ], [ 'CG2', 'HG21' ] ]
    },
    'LEU' => {
	'HB2' => [ [ 'CA', 'CB' ], [ 'CB', 'CG' ] ],
	'HB3' => [ [ 'CA', 'CB' ], [ 'CB', 'CG' ] ],
	'CG'  => [ [ 'CA', 'CB' ], [ 'CB', 'CG' ] ],
	'HG'  => [ [ 'CA', 'CB' ], [ 'CB', 'CG' ], [ 'CG', 'CD1' ] ],
	'CD1' => [ [ 'CA', 'CB' ], [ 'CB', 'CG' ], [ 'CG', 'CD1' ] ],
	'HD11'=> [ [ 'CA', 'CB' ], [ 'CB', 'CG' ], [ 'CG', 'CD1' ], [ 'CD1', 'HD11' ] ],
	'HD12'=> [ [ 'CA', 'CB' ], [ 'CB', 'CG' ], [ 'CG', 'CD1' ], [ 'CD1', 'HD11' ] ],
	'HD13'=> [ [ 'CA', 'CB' ], [ 'CB', 'CG' ], [ 'CG', 'CD1' ], [ 'CD1', 'HD11' ] ],
	'CD2' => [ [ 'CA', 'CB' ], [ 'CB', 'CG' ], [ 'CG', 'CD1' ] ],
	'HD21'=> [ [ 'CA', 'CB' ], [ 'CB', 'CG' ], [ 'CG', 'CD2' ], [ 'CD2', 'HD21' ] ],
	'HD22'=> [ [ 'CA', 'CB' ], [ 'CB', 'CG' ], [ 'CG', 'CD2' ], [ 'CD2', 'HD21' ] ],
	'HD23'=> [ [ 'CA', 'CB' ], [ 'CB', 'CG' ], [ 'CG', 'CD2' ], [ 'CD2', 'HD21' ] ]
    },
    'MET' => {
	'HB2' => [ [ 'CA', 'CB' ], [ 'CB', 'CG' ] ],
	'HB3' => [ [ 'CA', 'CB' ], [ 'CB', 'CG' ] ],
	'CG'  => [ [ 'CA', 'CB' ], [ 'CB', 'CG' ] ],
	'HG2' => [ [ 'CA', 'CB' ], [ 'CB', 'CG' ], [ 'CG', 'CD' ] ],
	'HG3' => [ [ 'CA', 'CB' ], [ 'CB', 'CG' ], [ 'CG', 'CD' ] ],
	'SD'  => [ [ 'CA', 'CB' ], [ 'CB', 'CG' ], [ 'CG', 'SD' ] ],
	'CE'  => [ [ 'CA', 'CB' ], [ 'CB', 'CG' ], [ 'CG', 'SD' ], [ 'SD', 'CE' ] ],
	'HE1' => [ [ 'CA', 'CB' ], [ 'CB', 'CG' ], [ 'CG', 'SD' ], [ 'SD', 'CE' ], [ 'CE', 'HE1' ] ],
	'HE2' => [ [ 'CA', 'CB' ], [ 'CB', 'CG' ], [ 'CG', 'SD' ], [ 'SD', 'CE' ], [ 'CE', 'HE1' ] ],
	'HE3' => [ [ 'CA', 'CB' ], [ 'CB', 'CG' ], [ 'CG', 'SD' ], [ 'SD', 'CE' ], [ 'CE', 'HE1' ] ]
    },
    'PHE' => {
	'HB2' => [ [ 'CA', 'CB' ], [ 'CB', 'CG' ] ],
	'HB3' => [ [ 'CA', 'CB' ], [ 'CB', 'CG' ] ],
	'CG'  => [ [ 'CA', 'CB' ], [ 'CB', 'CG' ] ],
	'HD1' => [ [ 'CA', 'CB' ], [ 'CB', 'CG' ], [ 'CG', 'CD1' ] ],
	'CD1' => [ [ 'CA', 'CB' ], [ 'CB', 'CG' ], [ 'CG', 'CD1' ] ],
	'HD2' => [ [ 'CA', 'CB' ], [ 'CB', 'CG' ], [ 'CG', 'CD1' ] ],
	'CD2' => [ [ 'CA', 'CB' ], [ 'CB', 'CG' ], [ 'CG', 'CD1' ] ],
	'HE1' => [ [ 'CA', 'CB' ], [ 'CB', 'CG' ], [ 'CG', 'CD1' ] ],
	'CE1' => [ [ 'CA', 'CB' ], [ 'CB', 'CG' ], [ 'CG', 'CD1' ] ],
	'HE2' => [ [ 'CA', 'CB' ], [ 'CB', 'CG' ], [ 'CG', 'CD1' ] ],
	'CE2' => [ [ 'CA', 'CB' ], [ 'CB', 'CG' ], [ 'CG', 'CD1' ] ],
	'HZ'  => [ [ 'CA', 'CB' ], [ 'CB', 'CG' ], [ 'CG', 'CD1' ] ],
	'CZ'  => [ [ 'CA', 'CB' ], [ 'CB', 'CG' ], [ 'CG', 'CD1' ] ]
    },
    'TYR' => {
	'HB2' => [ [ 'CA', 'CB' ], [ 'CB', 'CG' ] ],
	'HB3' => [ [ 'CA', 'CB' ], [ 'CB', 'CG' ] ],
	'CG'  => [ [ 'CA', 'CB' ], [ 'CB', 'CG' ] ],
	'HD1' => [ [ 'CA', 'CB' ], [ 'CB', 'CG' ], [ 'CG', 'CD1' ] ],
	'CD1' => [ [ 'CA', 'CB' ], [ 'CB', 'CG' ], [ 'CG', 'CD1' ] ],
	'HD2' => [ [ 'CA', 'CB' ], [ 'CB', 'CG' ], [ 'CG', 'CD1' ] ],
	'CD2' => [ [ 'CA', 'CB' ], [ 'CB', 'CG' ], [ 'CG', 'CD1' ] ],
	'HE1' => [ [ 'CA', 'CB' ], [ 'CB', 'CG' ], [ 'CG', 'CD1' ] ],
	'CE1' => [ [ 'CA', 'CB' ], [ 'CB', 'CG' ], [ 'CG', 'CD1' ] ],
	'HE2' => [ [ 'CA', 'CB' ], [ 'CB', 'CG' ], [ 'CG', 'CD1' ] ],
	'CE2' => [ [ 'CA', 'CB' ], [ 'CB', 'CG' ], [ 'CG', 'CD1' ] ],
	'CZ'  => [ [ 'CA', 'CB' ], [ 'CB', 'CG' ], [ 'CD1' ] ],
	'OH'  => [ [ 'CA', 'CB' ], [ 'CB', 'CG' ], [ 'CG', 'CD1' ] ],
	'HH'  => [ [ 'CA', 'CB' ], [ 'CB', 'CG' ], [ 'CG', 'CZ' ], [ 'CZ', 'OH' ], [ 'OH', 'HH' ] ]
    },
    'TRP' => {
	'HB2' => [ [ 'CA', 'CB' ], [ 'CB', 'CG' ] ],
	'HB3' => [ [ 'CA', 'CB' ], [ 'CB', 'CG' ] ],
	'CG'  => [ [ 'CA', 'CB' ], [ 'CB', 'CG' ] ],
	'HD1' => [ [ 'CA', 'CB' ], [ 'CB', 'CG' ], [ 'CG', 'CD1' ] ],
	'CD1' => [ [ 'CA', 'CB' ], [ 'CB', 'CG' ], [ 'CG', 'CD1' ] ],
	'CD2' => [ [ 'CA', 'CB' ], [ 'CB', 'CG' ], [ 'CG', 'CD1' ] ],
	'HE1' => [ [ 'CA', 'CB' ], [ 'CB', 'CG' ], [ 'CG', 'CD1' ] ],
	'NE1' => [ [ 'CA', 'CB' ], [ 'CB', 'CG' ], [ 'CG', 'CD1' ] ],
	'CE2' => [ [ 'CA', 'CB' ], [ 'CB', 'CG' ], [ 'CG', 'CD1' ] ],
	'HE3' => [ [ 'CA', 'CB' ], [ 'CB', 'CG' ], [ 'CG', 'CD1' ] ],
	'CE3' => [ [ 'CA', 'CB' ], [ 'CB', 'CG' ], [ 'CG', 'CD1' ] ],
	'HZ2' => [ [ 'CA', 'CB' ], [ 'CB', 'CG' ], [ 'CG', 'CD1' ] ],
	'CZ2' => [ [ 'CA', 'CB' ], [ 'CB', 'CG' ], [ 'CG', 'CD1' ] ],
	'HZ3' => [ [ 'CA', 'CB' ], [ 'CB', 'CG' ], [ 'CG', 'CD1' ] ],
	'CZ3' => [ [ 'CA', 'CB' ], [ 'CB', 'CG' ], [ 'CG', 'CD1' ] ],
	'HH2' => [ [ 'CA', 'CB' ], [ 'CB', 'CG' ], [ 'CG', 'CD1' ] ],
	'CH2' => [ [ 'CA', 'CB' ], [ 'CB', 'CG' ], [ 'CG', 'CD1' ] ]
    },
    'THR' => {
	'OG1' => [ [ 'CA', 'CB' ], [ 'CB', 'OG1' ] ],
	'HG1' => [ [ 'CA', 'CB' ], [ 'CB', 'OG1' ], [ 'OG1', 'HG1' ] ],
	'CG2' => [ [ 'CA', 'CB' ], [ 'CB', 'OG1' ] ],
	'HG21'=> [ [ 'CA', 'CB' ], [ 'CB', 'CG2' ], [ 'CG2', 'HG21' ] ],
	'HG22'=> [ [ 'CA', 'CB' ], [ 'CB', 'CG2' ], [ 'CG2', 'HG21' ] ],
	'HG23'=> [ [ 'CA', 'CB' ], [ 'CB', 'CG2' ], [ 'CG2', 'HG21' ] ]
    },
    'ASN' => {
	'HB2' => [ [ 'CA', 'CB' ], [ 'CB', 'CG' ] ],
	'HB3' => [ [ 'CA', 'CB' ], [ 'CB', 'CG' ] ],
	'CG'  => [ [ 'CA', 'CB' ], [ 'CB', 'CG' ] ],
	'HD21'=> [ [ 'CA', 'CB' ], [ 'CB', 'CG' ], [ 'CG', 'OD1' ] ],
	'ND2' => [ [ 'CA', 'CB' ], [ 'CB', 'CG' ], [ 'CG', 'OD1' ] ],
	'OD1' => [ [ 'CA', 'CB' ], [ 'CB', 'CG' ], [ 'CG', 'OD1' ] ]
    },
    'GLN' => {
	'HB2' => [ [ 'CA', 'CB' ], [ 'CB', 'CG' ] ],
	'HB3' => [ [ 'CA', 'CB' ], [ 'CB', 'CG' ] ],
	'CG'  => [ [ 'CA', 'CB' ], [ 'CB', 'CG' ] ],
	'HG2' => [ [ 'CA', 'CB' ], [ 'CB', 'CG' ], [ 'CG', 'CD' ] ],
	'HG3' => [ [ 'CA', 'CB' ], [ 'CB', 'CG' ], [ 'CG', 'CD' ] ],
	'CD'  => [ [ 'CA', 'CB' ], [ 'CB', 'CG' ], [ 'CG', 'CD' ] ],
	'HE21'=> [ [ 'CA', 'CB' ], [ 'CB', 'CG' ], [ 'CG', 'CD' ], [ 'CD', 'OE1' ] ],
	'NE2' => [ [ 'CA', 'CB' ], [ 'CB', 'CG' ], [ 'CG', 'CD' ], [ 'CD', 'OE1' ] ],
	'OE1' => [ [ 'CA', 'CB' ], [ 'CB', 'CG' ], [ 'CG', 'CD' ], [ 'CD', 'OE1' ] ]
    },
);

our %BOND_TYPES = (
    'single' => {
	'H' => {
	    'H' => {
		'min_length' => 2 * ( $ATOMS{'H'}{'covalent_radius'}{'length'}->[0]
				    - $ATOMS{'H'}{'covalent_radius'}{'error'}->[0] ),
		'max_length' => 2 * ( $ATOMS{'H'}{'covalent_radius'}{'length'}->[0]
				    + $ATOMS{'H'}{'covalent_radius'}{'error'}->[0] )
	    },
	    'C' => {
		'min_length' => $ATOMS{'H'}{'covalent_radius'}{'length'}->[0]
		              + $ATOMS{'C'}{'covalent_radius'}{'length'}->[0]
		              - $ATOMS{'H'}{'covalent_radius'}{'error'}->[0]
		              - $ATOMS{'C'}{'covalent_radius'}{'error'}->[0],
		'max_length' => $ATOMS{'H'}{'covalent_radius'}{'length'}->[0]
		              + $ATOMS{'C'}{'covalent_radius'}{'length'}->[0]
		              + $ATOMS{'H'}{'covalent_radius'}{'error'}->[0]
		              + $ATOMS{'C'}{'covalent_radius'}{'error'}->[0]
	    },
	    'N' => {
		'min_length' => $ATOMS{'H'}{'covalent_radius'}{'length'}->[0]
		              + $ATOMS{'N'}{'covalent_radius'}{'length'}->[0]
		              - $ATOMS{'H'}{'covalent_radius'}{'error'}->[0]
		              - $ATOMS{'N'}{'covalent_radius'}{'error'}->[0],
		'max_length' => $ATOMS{'H'}{'covalent_radius'}{'length'}->[0]
		              + $ATOMS{'N'}{'covalent_radius'}{'length'}->[0]
		              + $ATOMS{'H'}{'covalent_radius'}{'error'}->[0]
		              + $ATOMS{'N'}{'covalent_radius'}{'error'}->[0]
	    },
	    'O' => {
		'min_length' => $ATOMS{'H'}{'covalent_radius'}{'length'}->[0]
		              + $ATOMS{'O'}{'covalent_radius'}{'length'}->[0]
		              - $ATOMS{'H'}{'covalent_radius'}{'error'}->[0]
		              - $ATOMS{'O'}{'covalent_radius'}{'error'}->[0],
		'max_length' => $ATOMS{'H'}{'covalent_radius'}{'length'}->[0]
		              + $ATOMS{'O'}{'covalent_radius'}{'length'}->[0]
		              + $ATOMS{'H'}{'covalent_radius'}{'error'}->[0]
		              + $ATOMS{'O'}{'covalent_radius'}{'error'}->[0]
	    },
	    'S' => {
		'min_length' => $ATOMS{'H'}{'covalent_radius'}{'length'}->[0]
		              + $ATOMS{'S'}{'covalent_radius'}{'length'}->[0]
		              - $ATOMS{'H'}{'covalent_radius'}{'error'}->[0]
		              - $ATOMS{'S'}{'covalent_radius'}{'error'}->[0],
		'max_length' => $ATOMS{'H'}{'covalent_radius'}{'length'}->[0]
		              + $ATOMS{'S'}{'covalent_radius'}{'length'}->[0]
		              + $ATOMS{'H'}{'covalent_radius'}{'error'}->[0]
		              + $ATOMS{'S'}{'covalent_radius'}{'error'}->[0]
	    }
	},
	'C' => {
	    'C' => {
		'min_length' => 2 * ( $ATOMS{'C'}{'covalent_radius'}{'length'}->[0]
				    - $ATOMS{'C'}{'covalent_radius'}{'error'}->[0] ),
		'max_length' => 2 * ( $ATOMS{'C'}{'covalent_radius'}{'length'}->[0]
				    + $ATOMS{'C'}{'covalent_radius'}{'error'}->[0] )
	    },
	    'N' => {
		'min_length' => $ATOMS{'C'}{'covalent_radius'}{'length'}->[0]
		              + $ATOMS{'N'}{'covalent_radius'}{'length'}->[0]
		              - $ATOMS{'C'}{'covalent_radius'}{'error'}->[0]
		              - $ATOMS{'N'}{'covalent_radius'}{'error'}->[0],
		'max_length' => $ATOMS{'C'}{'covalent_radius'}{'length'}->[0]
		              + $ATOMS{'N'}{'covalent_radius'}{'length'}->[0]
		              + $ATOMS{'C'}{'covalent_radius'}{'error'}->[0]
		              + $ATOMS{'N'}{'covalent_radius'}{'error'}->[0]
	    },
	    'O' => {
		'min_length' => $ATOMS{'C'}{'covalent_radius'}{'length'}->[0]
		              + $ATOMS{'O'}{'covalent_radius'}{'length'}->[0]
		              - $ATOMS{'C'}{'covalent_radius'}{'error'}->[0]
		              - $ATOMS{'O'}{'covalent_radius'}{'error'}->[0],
		'max_length' => $ATOMS{'C'}{'covalent_radius'}{'length'}->[0]
		              + $ATOMS{'O'}{'covalent_radius'}{'length'}->[0]
		              + $ATOMS{'C'}{'covalent_radius'}{'error'}->[0]
		              + $ATOMS{'O'}{'covalent_radius'}{'error'}->[0]
	    },
	    'S' => {
		'min_length' => $ATOMS{'C'}{'covalent_radius'}{'length'}->[0]
		              + $ATOMS{'S'}{'covalent_radius'}{'length'}->[0]
		              - $ATOMS{'C'}{'covalent_radius'}{'error'}->[0]
		              - $ATOMS{'S'}{'covalent_radius'}{'error'}->[0],
		'max_length' => $ATOMS{'C'}{'covalent_radius'}{'length'}->[0]
		              + $ATOMS{'S'}{'covalent_radius'}{'length'}->[0]
		              + $ATOMS{'C'}{'covalent_radius'}{'error'}->[0]
		              + $ATOMS{'S'}{'covalent_radius'}{'error'}->[0]
	    }
	},
	'N' => {
	    'N' => {
		'min_length' => 2 * ( $ATOMS{'N'}{'covalent_radius'}{'length'}->[0]
				    - $ATOMS{'N'}{'covalent_radius'}{'error'}->[0] ),
		'max_length' => 2 * ( $ATOMS{'N'}{'covalent_radius'}{'length'}->[0]
				    + $ATOMS{'N'}{'covalent_radius'}{'error'}->[0] )
	    },
	    'O' => {
		'min_length' => $ATOMS{'N'}{'covalent_radius'}{'length'}->[0]
		              + $ATOMS{'O'}{'covalent_radius'}{'length'}->[0]
		              - $ATOMS{'N'}{'covalent_radius'}{'error'}->[0]
		              - $ATOMS{'O'}{'covalent_radius'}{'error'}->[0],
		'max_length' => $ATOMS{'N'}{'covalent_radius'}{'length'}->[0]
		              + $ATOMS{'O'}{'covalent_radius'}{'length'}->[0]
		              + $ATOMS{'N'}{'covalent_radius'}{'error'}->[0]
		              + $ATOMS{'O'}{'covalent_radius'}{'error'}->[0]
	    },
	    'S' => {
		'min_length' => $ATOMS{'N'}{'covalent_radius'}{'length'}->[0]
		              + $ATOMS{'S'}{'covalent_radius'}{'length'}->[0]
		              - $ATOMS{'N'}{'covalent_radius'}{'error'}->[0]
		              - $ATOMS{'S'}{'covalent_radius'}{'error'}->[0],
		'max_length' => $ATOMS{'N'}{'covalent_radius'}{'length'}->[0]
		              + $ATOMS{'S'}{'covalent_radius'}{'length'}->[0]
		              + $ATOMS{'N'}{'covalent_radius'}{'error'}->[0]
		              + $ATOMS{'S'}{'covalent_radius'}{'error'}->[0]
	    }
	},
	'O' => {
	    'O' => {
		'min_length' => 2 * ( $ATOMS{'O'}{'covalent_radius'}{'length'}->[0]
				    - $ATOMS{'O'}{'covalent_radius'}{'error'}->[0] ),
		'max_length' => 2 * ( $ATOMS{'O'}{'covalent_radius'}{'length'}->[0]
				    + $ATOMS{'O'}{'covalent_radius'}{'error'}->[0] )
	    },
	    'S' => {
		'min_length' => $ATOMS{'O'}{'covalent_radius'}{'length'}->[0]
		              + $ATOMS{'S'}{'covalent_radius'}{'length'}->[0]
		              - $ATOMS{'O'}{'covalent_radius'}{'error'}->[0]
		              - $ATOMS{'S'}{'covalent_radius'}{'error'}->[0],
		'max_length' => $ATOMS{'O'}{'covalent_radius'}{'length'}->[0]
		              + $ATOMS{'S'}{'covalent_radius'}{'length'}->[0]
		              + $ATOMS{'O'}{'covalent_radius'}{'error'}->[0]
		              + $ATOMS{'S'}{'covalent_radius'}{'error'}->[0]
	    }
	},
	'S' => {
	    'S' => {
		'min_length' => 2 * ( $ATOMS{'S'}{'covalent_radius'}{'length'}->[0]
				    - $ATOMS{'S'}{'covalent_radius'}{'error'}->[0] ),
		'max_length' => 2 * ( $ATOMS{'S'}{'covalent_radius'}{'length'}->[0]
				    + $ATOMS{'S'}{'covalent_radius'}{'error'}->[0] )
	    }
	}
    },
    'double' => {
	'C' => {
	    'C' => {
		'min_length' => 2 * ( $ATOMS{'C'}{'covalent_radius'}{'length'}->[1]
				    - $ATOMS{'C'}{'covalent_radius'}{'error'}->[1] ),
		'max_length' => 2 * ( $ATOMS{'C'}{'covalent_radius'}{'length'}->[1]
				    + $ATOMS{'C'}{'covalent_radius'}{'error'}->[1] )
	    },
	    'N' => {
		'min_length' => $ATOMS{'C'}{'covalent_radius'}{'length'}->[1]
		              + $ATOMS{'N'}{'covalent_radius'}{'length'}->[1]
		              - $ATOMS{'C'}{'covalent_radius'}{'error'}->[1]
		              - $ATOMS{'N'}{'covalent_radius'}{'error'}->[1],
		'max_length' => $ATOMS{'C'}{'covalent_radius'}{'length'}->[1]
		              + $ATOMS{'N'}{'covalent_radius'}{'length'}->[1]
		              + $ATOMS{'C'}{'covalent_radius'}{'error'}->[1]
		              + $ATOMS{'N'}{'covalent_radius'}{'error'}->[1]
	    },
	    'O' => {
		'min_length' => $ATOMS{'C'}{'covalent_radius'}{'length'}->[1]
		              + $ATOMS{'O'}{'covalent_radius'}{'length'}->[1]
		              - $ATOMS{'C'}{'covalent_radius'}{'error'}->[1]
		              - $ATOMS{'O'}{'covalent_radius'}{'error'}->[1],
		'max_length' => $ATOMS{'C'}{'covalent_radius'}{'length'}->[1]
		              + $ATOMS{'O'}{'covalent_radius'}{'length'}->[1]
		              + $ATOMS{'C'}{'covalent_radius'}{'error'}->[1]
		              + $ATOMS{'O'}{'covalent_radius'}{'error'}->[1]
	    }
	},
	'N' => {
	    'N' => {
		'min_length' => 2 * ( $ATOMS{'N'}{'covalent_radius'}{'length'}->[1]
				    - $ATOMS{'N'}{'covalent_radius'}{'error'}->[1] ),
		'max_length' => 2 * ( $ATOMS{'N'}{'covalent_radius'}{'length'}->[1]
				    + $ATOMS{'N'}{'covalent_radius'}{'error'}->[1] )
	    },
	    'O' => {
		'min_length' => $ATOMS{'N'}{'covalent_radius'}{'length'}->[1]
		              + $ATOMS{'O'}{'covalent_radius'}{'length'}->[1]
		              - $ATOMS{'N'}{'covalent_radius'}{'error'}->[1]
		              - $ATOMS{'O'}{'covalent_radius'}{'error'}->[1],
		'max_length' => $ATOMS{'N'}{'covalent_radius'}{'length'}->[1]
		              + $ATOMS{'O'}{'covalent_radius'}{'length'}->[1]
		              + $ATOMS{'N'}{'covalent_radius'}{'error'}->[1]
		              + $ATOMS{'O'}{'covalent_radius'}{'error'}->[1]
	    }
	},
	'O' => {
	    'O' => {
		'min_length' => 2 * ( $ATOMS{'O'}{'covalent_radius'}{'length'}->[1]
				    - $ATOMS{'O'}{'covalent_radius'}{'error'}->[1] ),
		'max_length' => 2 * ( $ATOMS{'O'}{'covalent_radius'}{'length'}->[1]
				    + $ATOMS{'O'}{'covalent_radius'}{'error'}->[1] )
	    }
	}
    },
    'triple' => {
	'C' => {
	    'C' => {
		'min_length' => 2 * ( $ATOMS{'C'}{'covalent_radius'}{'length'}->[2]
				    - $ATOMS{'C'}{'covalent_radius'}{'error'}->[2] ),
		'max_length' => 2 * ( $ATOMS{'C'}{'covalent_radius'}{'length'}->[2]
				    + $ATOMS{'C'}{'covalent_radius'}{'error'}->[2] )
	    }
	}
    }
);

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
