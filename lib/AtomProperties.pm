package AtomProperties;

use strict;
use warnings;

use Exporter qw( import );
our @EXPORT_OK = qw( %ATOMS
                     %HYDROGEN_NAMES
                     @MAINCHAIN_NAMES
                     @MAINCHAIN_GRAMMAR
                     @SIDECHAIN_GRAMMAR
                     sort_atom_names );

use Version qw( $VERSION );

our $VERSION = $VERSION;

# ----------------------------- Atom properties ------------------------------- #

our %ATOMS = (
    q{.} => { # Point in space.
             'covalent_radius' => {
                                    'length'    => [ 0 ],
                                    'error'     => [ 0 ]
                                  },
             'lone_pairs' => 0,
             'vdw_radius' => 0,
             'valence' => 0,
             'partial_charge' => 0
           },
    'X' => { # Dummy atom.
             'covalent_radius' => {
                                    'length'    => [ 1, 0.5, 0.25 ],
                                    'error'     => [ 0.05, 0.05, 0.05 ]
                                  },
             'lone_pairs' => 0,
             'vdw_radius' => 0,
             'valence' => 0,
             'partial_charge' => 0
           },
    'H' => {
             'covalent_radius' => {
                                    'length'    => [ 0.31 ],
                                    'error'     => [ 0.05 ]
                                  },
             'lone_pairs' => 0,
             'vdw_radius' => 1.2,
             'valence' => 1,
             'partial_charge' => 0.27
           },
    'C' => {
             'covalent_radius' => {
                                    'length'    => [ 0.76, 0.70, 0.64 ],
                                    'error'     => [ 0.03, 0.03, 0.03 ]
                                  },
             'lone_pairs' => 0,
             'vdw_radius' => 1.77,
             'valence' => 4,
             'partial_charge' => 0.12,
           },
    'N' => {
             'covalent_radius' => {
                                    'length'    => [ 0.71, 0.60 ],
                                    'error'     => [ 0.02, 0.06 ]
                                  },
             'lone_pairs' => 1,
             'vdw_radius' => 1.66,
             'valence' => 3,
             'partial_charge' => -0.42
           },
    'O' => {
             'covalent_radius' => {
                                    'length'    => [ 0.66, 0.57 ],
                                    'error'     => [ 0.02, 0.07 ]
                                  },
             'lone_pairs' => 2,
             'vdw_radius' => 1.5,
             'valence' => 2,
             'partial_charge' => -0.57,
           },
    'S' => {
             'covalent_radius' => {
                                    'length'    => [ 1.05 ],
                                    'error'     => [ 0.03 ]
                                  },
             'lone_pairs' => 2,
             'vdw_radius' => 1.89,
             'valence' => 2,
             'partial_charge' => 0
           },
    'P' => {
             'covalent_radius' => {
                                    'length'    => [ 1.07 ],
                                    'error'     => [ 0.03 ]
                                  },
             'lone_pairs' => 0,
             'vdw_radius' => 1.9,
             'valence' => 5,
             'partial_charge' => 0
           },
    );

# -------------------------------- Atom names --------------------------------- #

our %HYDROGEN_NAMES = (
    'XAA' => { # Dummy side-chain.
        'CA'  => [ 'HA' ],
        'CB'  => [ 'HB2', 'HB3' ],
        'OG'  => [ 'HG' ]
    },
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
    },
    'GLU' => {
        'N'   => [ 'H', 'H2' ],
        'CA'  => [ 'HA' ],
        'OXT' => [ 'HXT' ],
        'CB'  => [ 'HB2', 'HB3' ],
        'CG'  => [ 'HG2', 'HG3' ],
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
        'N'   => [ 'H', 'H2' ],
        'CA'  => [ 'HA' ],
        'OXT' => [ 'HXT' ],
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
    },
);

# TODO: should include proline atoms.   v should not be here
our @MAINCHAIN_NAMES = qw( N CA C O OXT CB H H2 HA HA2 HA3 HXT XA XA2 );

our @MAINCHAIN_GRAMMAR = qw( N CA C O OXT H H2 HA HA2 HA3 HXT XA XA2 );

our @SIDECHAIN_GRAMMAR = qw( CB CD CD1 CD2 CE CE1 CE2 CE3 CG CG1 CG2 CH2 CZ CZ2
                             CZ3 HB HB1 HB2 HB3 HD1 HD11 HD12 HD13 HD2 HD21
                             HD22 HD23 HD3 HE HE1 HE2 HE21 HE22 HE3 HG HG1 HG11
                             HG12 HG13 HG2 HG21 HG22 HG23 HG3 HH HH11 HH12 HH2
                             HH21 HH22 HZ HZ1 HZ2 HZ3 ND1 ND2 NE NE1 NE2 NH1 NH2
                             NZ OD1 OD2 OE1 OE2 OG OG1 OH SD SG );

#
# Sorts atom names by their hierarchical rules.
# Input:
#     $atom_names - list of atom names;
#     $options{sort_type} - sorts by atom type, greek letter, number or their
#     combinations.
# Output:
#     @sorted_names - sorted list of atom names.
#

sub sort_atom_names
{
    my ( $atom_names, $options ) = @_;
    my ( $sort_type ) = ( $options->{'sort_type'} );

    $sort_type //= 'tgn';

    # First priority is decided by atom type: S > P > O > N > C > H.
    # Second priority - by greek letter: A > B > G > D > E > Z > H.
    # TODO: look for more greek letters that are in PDBx.
    # Third priority - by numeration: 1 > 2 > 3 and etc.
    # This priority is achieved by assinging first, second and third priorities
    # to numbers. Then iteratively is sorted by priorities.
    my %atom_type_priority =
        ( 'H' => 1, 'X' => 2, 'C' => 3, 'N' => 4, 'O' => 5, 'P' => 6, 'S' => 7,);
    my %greek_letter_priority =
        ( 'H' => 1, 'Z' => 2, 'E' => 3, 'D' => 4, 'G' => 5, 'B' => 6, 'A' => 7,
          q() => 8, );

    # Decomposes each atom name by its components.
    my %atom_names;
    for my $atom_name ( @{ $atom_names } ) {
        my ( $atom_type ) =
            $atom_name =~ /(^\p{IsAlphabetic})\p{IsAlphabetic}?\d?/smx;
        my ( $greek_letter ) =
            $atom_name =~ /^\p{IsAlphabetic}(\p{IsAlphabetic}?)\d?/smx;
        my ( $number ) =
            $atom_name =~ /^\p{IsAlphabetic}\p{IsAlphabetic}?(\d?)/smx;
        $atom_names{$atom_name}{'type'} =
            $atom_type_priority{$atom_type};
        $atom_names{$atom_name}{'greek_letter'} =
            $greek_letter_priority{$greek_letter};
        $atom_names{$atom_name}{'number'} = $number;
    }

    # Sorts by rules of described in %atom_names.
    my @sorted_names;
    if( $sort_type eq 'tgn') { # By type, then greek letter and then number.
        @sorted_names =
          sort {
            $atom_names{$b}{'type'} <=> $atom_names{$a}{'type'} ||
            $atom_names{$b}{'greek_letter'} <=> $atom_names{$a}{'greek_letter'}||
            $atom_names{$a}{'number'} cmp $atom_names{$b}{'number'} }
              @{ $atom_names };
    } elsif( $sort_type eq 'gn' ) { # By greek letter and then number.
        @sorted_names =
          sort {
            $atom_names{$b}{'greek_letter'} <=> $atom_names{$a}{'greek_letter'}||
            $atom_names{$a}{'number'} cmp $atom_names{$b}{'number'} }
              @{ $atom_names };
    }

    return \@sorted_names;
}

1;
