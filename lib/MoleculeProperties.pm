package MoleculeProperties;

use strict;
use warnings;

use Exporter qw( import );
our @EXPORT_OK = qw( %BOND_TYPES
                     %HYDROGEN_NAMES );

use AtomProperties qw( %ATOMS );

# ------------------------------- Molecule properties ------------------------- #

#
# Stores molecule properties, such as hybridization angles, hydrogen names for
# specific residues.
#

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
