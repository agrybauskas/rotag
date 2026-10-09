#!/usr/bin/perl

use strict;
use warnings;

use LinearAlgebra qw( mult_matrix_product );
use Symbolic;

my $matrices = [
  [
    [
      '0.537003926704602',
      '-0.678207399996524',
      '-0.501658753829526',
      '34.297'
    ],
    [
      '0.409911012029786',
      '-0.309965098545599',
      '0.857842992569347',
      '2.104'
    ],
    [
      '-0.737292170662726',
      '-0.666300502981629',
      '0.111552206638299',
      '42.669'
    ],
    [
      0,
      0,
      0,
      1
    ]
  ],
  bless( {
    'is_evaluated' => undef,
    'matrix' => sub {
        package AlterMolecule;
        use warnings;
        use strict;
        my($svar) = @_;
        return [[cos $svar, -sin($svar), 0, 0], [sin $svar, cos $svar, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]];
    },
    'symbols' => [
      'chi1'
    ]
  }, 'Symbolic' ),
  [
    [
      '-0.947029553920039',
      '0.123452632652295',
      '0.296470017865601',
      '-3.5527136788005e-15'
    ],
    [
      '-0.321146421437345',
      '-0.364049803537254',
      '-0.874261240443881',
      '7.105427357601e-15'
    ],
    [
      '-1.80411241501588e-16',
      '-0.923161517848152',
      '0.384412294241867',
      '1.53291454425874'
    ],
    [
      0,
      0,
      0,
      1
    ]
  ],
  bless( {
    'is_evaluated' => undef,
    'matrix' => sub {
        package AlterMolecule;
        use warnings;
        use strict;
        my($svar) = @_;
        return [[cos $svar, -sin($svar), 0, 0], [sin $svar, cos $svar, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]];
    },
    'symbols' => [
      'chi2'
    ]
  }, 'Symbolic' ),
  [
    [
      '-0.947029553920039',
      '-0.321146421437345',
      '-1.80411241501588e-16',
      '-3.5527136788005e-15'
    ],
    [
      '0.123452632652295',
      '-0.364049803537254',
      '-0.923161517848152',
      '1.41512771740941'
    ],
    [
      '0.296470017865601',
      '-0.874261240443881',
      '0.384412294241867',
      '-0.589271196835234'
    ],
    [
      0,
      0,
      0,
      1
    ]
  ],
  bless( {
    'is_evaluated' => undef,
    'matrix' => sub {
        package AlterMolecule;
        use warnings;
        use strict;
        my($svar) = @_;
        return [[cos $svar, 0, sin $svar, 0], [0, 1, 0, 0], [-sin($svar), 0, cos $svar, 0], [0, 0, 0, 1]];
    },
    'symbols' => [
      'eta'
    ]
  }, 'Symbolic' ),
  bless( {
    'is_evaluated' => undef,
    'matrix' => sub {
        package AlterMolecule;
        use warnings;
        use strict;
        my($svar) = @_;
        return [[1, 0, 0, 0], [0, cos $svar, -sin($svar), 0], [0, sin $svar, cos $svar, 0], [0, 0, 0, 1]];
    },
    'symbols' => [
      'N-CA-CB'
    ]
  }, 'Symbolic' ),
  [
    [
      '-0.947029553920039',
      '0.123452632652295',
      '0.296470017865601',
      '-3.5527136788005e-15'
    ],
    [
      '-0.321146421437345',
      '-0.364049803537254',
      '-0.874261240443881',
      '7.105427357601e-15'
    ],
    [
      '-1.80411241501588e-16',
      '-0.923161517848152',
      '0.384412294241867',
      '1.53291454425874'
    ],
    [
      0,
      0,
      0,
      1
    ]
  ],
  bless( {
    'is_evaluated' => undef,
    'matrix' => sub {
        package AlterMolecule;
        use warnings;
        use strict;
        my($svar) = @_;
        return [[cos $svar, 0, sin $svar, 0], [0, 1, 0, 0], [-sin($svar), 0, cos $svar, 0], [0, 0, 0, 1]];
    },
    'symbols' => [
      'eta'
    ]
  }, 'Symbolic' ),
  bless( {
    'is_evaluated' => undef,
    'matrix' => sub {
        package AlterMolecule;
        use warnings;
        use strict;
        my($svar) = @_;
        return [[1, 0, 0, 0], [0, cos $svar, -sin($svar), 0], [0, sin $svar, cos $svar, 0], [0, 0, 0, 1]];
    },
    'symbols' => [
      'CA-CB-CG'
    ]
  }, 'Symbolic' ),
  [
    [
      '0.420170108867224',
      '-0.426804926327117',
      '-0.8008087377629',
      '0'
    ],
    [
      '0.90744535902417',
      '0.197621455194589',
      '0.370794661278002',
      '0'
    ],
    [
      '2.77555756156289e-17',
      '-0.882487005745511',
      '0.470336777947804',
      '1.521565969651'
    ],
    [
      0,
      0,
      0,
      1
    ]
  ],
  bless( {
    'is_evaluated' => undef,
    'matrix' => sub {
        package AlterMolecule;
        use warnings;
        use strict;
        my($svar) = @_;
        return [[cos $svar, 0, sin $svar, 0], [0, 1, 0, 0], [-sin($svar), 0, cos $svar, 0], [0, 0, 0, 1]];
    },
    'symbols' => [
      'eta'
    ]
  }, 'Symbolic' ),
  bless( {
    'is_evaluated' => undef,
    'matrix' => sub {
        package AlterMolecule;
        use warnings;
        use strict;
        my($svar) = @_;
        return [[1, 0, 0, 0], [0, cos $svar, -sin($svar), 0], [0, sin $svar, cos $svar, 0], [0, 0, 0, 1]];
    },
    'symbols' => [
      'CB-CG-OD1'
    ]
  }, 'Symbolic' ),
  [
    [
      '-7.105427357601e-15'
    ],
    [
      '0'
    ],
    [
      '1.24913009730772'
    ],
    [
      1
    ]
  ]
];

my $matrix_product =
    mult_matrix_product( $matrices,
                         { 'CA-CB-CG' => '0.129027660387877',
                           'CB-CG-OD1' => '0',
                           'N-CA-CB' => '0',
                           'chi1' => '0',
                           'chi2' => '0.433632771010695',
                           'eta' => 0 } );

for my $matrix_id ( 0..$#{ $matrix_product } ) {
    for my $row ( @{ $matrix_product->[$matrix_id] } ) {
        print( join( " ", @{ $row } ), "\n" );
    }
    print( "\n" ) if $matrix_id != $#{ $matrix_product };
}
