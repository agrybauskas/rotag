package MoleculeProperties;

use strict;
use warnings;

use Exporter qw( import );
our @EXPORT_OK = qw( %HYBRIDIZATION
                     %ROTATABLE_BONDS );

# ------------------------------- Molecule properties ------------------------- #

#
# Stores molecule properties, such as rotatable bonds of certain sidechains.
#

our %ROTATABLE_BONDS = (
    'PHE' => {
	'CD2' => [ [ 'CA', 'CB' ], [ 'CB', 'CG' ], [ 'CG', 'CD1' ] ],
	'CE2' => [ [ 'CA', 'CB' ], [ 'CB', 'CG' ], [ 'CG', 'CD1' ] ],
	'CD1' => [ [ 'CA', 'CB' ], [ 'CB', 'CG' ], [ 'CG', 'CD1' ] ],
	'CZ'  => [ [ 'CA', 'CB' ], [ 'CB', 'CG' ], [ 'CG', 'CD1' ] ],
	'CE1' => [ [ 'CA', 'CB' ], [ 'CB', 'CG' ], [ 'CG', 'CD1' ] ],
	'CG'  => [ [ 'CA', 'CB' ], [ 'CB', 'CG' ] ]
    },
    'TRP' => {
	'CE3' => [ [ 'CA', 'CB' ], [ 'CB', 'CG' ], [ 'CD1' ] ],
	'CH2' => [ [ 'CA', 'CB' ], [ 'CB', 'CG' ], [ 'CD1' ] ],
	'CE2' => [ [ 'CA', 'CB' ], [ 'CB', 'CG' ], [ 'CD1' ] ],
	'CD2' => [ [ 'CA', 'CB' ], [ 'CB', 'CG' ], [ 'CD1' ] ],
	'CD1' => [ [ 'CA', 'CB' ], [ 'CB', 'CG' ], [ 'CD1' ] ],
	'CZ2' => [ [ 'CA', 'CB' ], [ 'CB', 'CG' ], [ 'CD1' ] ],
	'CZ3' => [ [ 'CA', 'CB' ], [ 'CB', 'CG' ], [ 'CD1' ] ],
	'NE1' => [ [ 'CA', 'CB' ], [ 'CB', 'CG' ], [ 'CD1' ] ],
	'CG'  => [ [ 'CA', 'CB' ], [ 'CB', 'CG' ] ]
    },
    'TYR' => {
	'CD1' => [ [ 'CA', 'CB' ], [ 'CB', 'CG' ], [ 'CD1' ] ],
	'CD2' => [ [ 'CA', 'CB' ], [ 'CB', 'CG' ], [ 'CD1' ] ],
	'CE2' => [ [ 'CA', 'CB' ], [ 'CB', 'CG' ], [ 'CD1' ] ],
	'CZ'  => [ [ 'CA', 'CB' ], [ 'CB', 'CG' ], [ 'CD1' ] ],
	'CG'  => [ [ 'CA', 'CB' ], [ 'CB', 'CG' ] ],
	'CE1' => [ [ 'CA', 'CB' ], [ 'CB', 'CG' ], [ 'CD1' ] ],
	'OH'  => [ [ 'CA', 'CB' ], [ 'CB', 'CG' ], [ 'CD1' ] ]
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
    'SER' => {
	'OG'  => [ [ 'CA', 'CB' ], [ 'CB', 'OG' ] ]
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
    'GLY' => { 'N'  => 'sp3',
	       'CA' => 'sp3',
	       'C'  => 'sp2',
	       'O'  => 'sp2',
	       'OXT'=> 'sp3'
    } );

1;
