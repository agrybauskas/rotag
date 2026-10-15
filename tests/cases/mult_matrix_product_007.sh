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
        return [[1, 0, 0, 0], [0, cos $svar, -sin($svar), 0], [0, sin $svar, cos $svar, 0], [0, 0, 0, 1]];
    },
    'symbols' => [
      'N-CA-CB'
    ]
  }, 'Symbolic' ),
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
  bless( {
    'is_evaluated' => undef,
    'matrix' => sub {
        package AlterMolecule;
        use warnings;
        use strict;
        my($svar) = @_;
        return [[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, $svar], [0, 0, 0, 1]];
    },
    'symbols' => [
      'CA-CB'
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
        return [[1, 0, 0, 0], [0, cos $svar, -sin($svar), 0], [0, sin $svar, cos $svar, 0], [0, 0, 0, 1]];
    },
    'symbols' => [
      'CA-CB-CG'
    ]
  }, 'Symbolic' ),
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
  bless( {
    'is_evaluated' => undef,
    'matrix' => sub {
        package AlterMolecule;
        use warnings;
        use strict;
        my($svar) = @_;
        return [[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, $svar], [0, 0, 0, 1]];
    },
    'symbols' => [
      'CB-CG'
    ]
  }, 'Symbolic' ),
  [
    [
      '-0.424479186506929',
      '0.422519568475836',
      '-0.8008087377629',
      '0'
    ],
    [
      '-0.905393416398365',
      '-0.206818956360942',
      '0.370794661278002',
      '0'
    ],
    [
      '-0.00895442711252087',
      '0.882441575145213',
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
        return [[cos $svar, -sin($svar), 0, 0], [sin $svar, cos $svar, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]];
    },
    'symbols' => [
      'CB-CG-OD1.H1{1}'
    ]
  }, 'Symbolic' ),
  bless( {
    'is_evaluated' => undef,
    'matrix' => sub {
        package AlterMolecule;
        use warnings;
        use strict;
        my($svar) = @_;
        return [[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, $svar], [0, 0, 0, 1]];
    },
    'symbols' => [
      'CG-OD1'
    ]
  }, 'Symbolic' ),
  [
    [
      '-0.573664118939681',
      '-0.0707741434565305',
      '0.816027266247369',
      '0'
    ],
    [
      '-0.819090641285297',
      '0.0495678800407138',
      '-0.571518631915977',
      '0'
    ],
    [
      '2.77555756156289e-17',
      '-0.996260029252536',
      '-0.0864057527814912',
      '1.24913009730772'
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
        return [[1, 0, 0, 0], [0, cos $svar, -sin($svar), 0], [0, sin $svar, cos $svar, 0], [0, 0, 0, 1]];
    },
    'symbols' => [
      'CG-OD1.H1{1}'
    ]
  }, 'Symbolic' ),
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
      'CG-OD1.H1{1}-O{1}'
    ]
  }, 'Symbolic' ),
  bless( {
    'is_evaluated' => undef,
    'matrix' => sub {
        package AlterMolecule;
        use warnings;
        use strict;
        my($svar) = @_;
        return [[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, $svar], [0, 0, 0, 1]];
    },
    'symbols' => [
      'OD1.H1{1}'
    ]
  }, 'Symbolic' ),
  [
    [
      '0.116552209485925',
      '-0.635424109637429',
      '0.763316306229226',
      '7.105427357601e-15'
    ],
    [
      '-0.99318456616278',
      '-0.0745682992487531',
      '0.0895767061417168',
      '-1.95399252334028e-14'
    ],
    [
      '-1.66533453693773e-16',
      '-0.768554337466538',
      '-0.63978451869467',
      '88.6577443317841'
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
        return [[1, 0, 0, 0], [0, cos $svar, -sin($svar), 0], [0, sin $svar, cos $svar, 0], [0, 0, 0, 1]];
    },
    'symbols' => [
      'OD1.H1{1}-O{1}'
    ]
  }, 'Symbolic' ),
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
      'OD1.H1{1}-O{1}.MG{2}'
    ]
  }, 'Symbolic' ),
  bless( {
    'is_evaluated' => undef,
    'matrix' => sub {
        package AlterMolecule;
        use warnings;
        use strict;
        my($svar) = @_;
        return [[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, $svar], [0, 0, 0, 1]];
    },
    'symbols' => [
      'H1{1}-O{1}'
    ]
  }, 'Symbolic' ),
  [
    [
      '0.999011011723179',
      '-0.0389096390586475',
      '-0.0215183280985074',
      '3.5527136788005e-15'
    ],
    [
      '0.0444634507863395',
      '0.874227195467821',
      '0.483476795977469',
      '0'
    ],
    [
      '-1.11022302462516e-16',
      '-0.483955422216545',
      '0.875092651841624',
      '0.970009278306144'
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
        return [[1, 0, 0, 0], [0, cos $svar, -sin($svar), 0], [0, sin $svar, cos $svar, 0], [0, 0, 0, 1]];
    },
    'symbols' => [
      'H1{1}-O{1}.MG{2}'
    ]
  }, 'Symbolic' ),
  bless( {
    'is_evaluated' => undef,
    'matrix' => sub {
        package AlterMolecule;
        use warnings;
        use strict;
        my($svar) = @_;
        return [[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, $svar], [0, 0, 0, 1]];
    },
    'symbols' => [
      'O{1}.MG{2}'
    ]
  }, 'Symbolic' ),
  [
    [
      '-1.21200845646577e-15'
    ],
    [
      '-1.29864511023586e-15'
    ],
    [
      '36.6089350568956'
    ],
    [
      1
    ]
  ]
];

my $matrix_product =
    mult_matrix_product( $matrices,
                         { 'CA-CB' => '0',
                           'CA-CB-CG' => '0',
                           'CB-CG' => '0',
                           'CB-CG-OD1' => '0',
                           'CB-CG-OD1.H1{1}' => '2.01738266762786',
                           'CG-OD1' => '0',
                           'CG-OD1.H1{1}' => '0.610112408113143',
                           'CG-OD1.H1{1}-O{1}' => '4.42103833711153',
                           'H1{1}-O{1}' => '0',
                           'H1{1}-O{1}.MG{2}' => '-0.542028454845521',
                           'N-CA-CB' => '0',
                           'OD1.H1{1}' => '-85.6577443317841',
                           'OD1.H1{1}-O{1}' => '0',
                           'OD1.H1{1}-O{1}.MG{2}' => '2.92258161385909',
                           'O{1}.MG{2}' => '-33.6089350568956',
                           'chi1' => '-4.70826062537455e-06',
                           'chi2' => '5.7183527037985e-06' } );

for my $matrix_id ( 0..$#{ $matrix_product } ) {
    for my $row ( @{ $matrix_product->[$matrix_id] } ) {
        print( join( " ", @{ $row } ), "\n" );
    }
    print( "\n" ) if $matrix_id != $#{ $matrix_product };
}
