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
        return [[cos $svar, -sin($svar), 0, 0], [sin $svar, cos $svar, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]];
    },
    'symbols' => [
      'CB-CG-OD1.H1{1}'
    ]
  }, 'Symbolic' ),
  [
    [
      '0.581945742910416',
      '0.0702675441776074',
      '-0.810186166596105',
      '-7.105427357601e-15'
    ],
    [
      '0.813227614083809',
      '-0.0502834599941585',
      '0.57976928285532',
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
        return [[cos $svar, -sin($svar), 0, 0], [sin $svar, cos $svar, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]];
    },
    'symbols' => [
      'CG-OD1.H1{1}-O{1}'
    ]
  }, 'Symbolic' ),
  [
    [
      '-0.0305919345607513',
      '0.409610081640015',
      '-0.911747615603513',
      '2.48960247807609'
    ],
    [
      '-0.96856085902541',
      '0.213162803827254',
      '0.128263328463215',
      '1.6932948677793'
    ],
    [
      '0.246888630568096',
      '0.887006877134199',
      '0.390211229993251',
      '0.350099974164053'
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
      '0.581945742910416',
      '0.0702675441776074',
      '-0.810186166596105',
      '-7.105427357601e-15'
    ],
    [
      '0.813227614083809',
      '-0.0502834599941585',
      '0.57976928285532',
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
      'CG-OD1.H1{1}'
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
      'OD1.H1{1}-O{1}'
    ]
  }, 'Symbolic' ),
  [
    [
      '0.958394139007897',
      '-0.163969046798359',
      '-0.233655357326446',
      '-1.39158565907965'
    ],
    [
      '-0.0980845391374959',
      '-0.957884291953252',
      '0.269883505259959',
      '66.1610035461018'
    ],
    [
      '-0.268067337617802',
      '-0.23573679161568',
      '-0.934113519643758',
      '58.5498976867266'
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
        return [[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, $svar], [0, 0, 0, 1]];
    },
    'symbols' => [
      'H1{1}-O{1}'
    ]
  }, 'Symbolic' ),
  [
    [
      '3.5527136788005e-15'
    ],
    [
      '0'
    ],
    [
      '0.970009278306144'
    ],
    [
      1
    ]
  ]
];

my $matrix_product =
    mult_matrix_product( $matrices,
    {
        'CA-CB' => '0',
        'CA-CB-CG' => '0',
        'CB-CG' => '0',
        'CB-CG-OD1' => '0',
        'CB-CG-OD1.H1{1}' => '2.01738266762786',
        'CG-OD1' => '0',
        'CG-OD1.H1{1}' => '0.610112408113143',
        'CG-OD1.H1{1}-O{1}' => '4.42103833711153',
        'H1{1}-O{1}' => '0',
        'N-CA-CB' => '0',
        'OD1.H1{1}' => '-85.6577443317841',
        'OD1.H1{1}-O{1}' => '0',
        'chi1' => '-4.70826062537455e-06',
        'chi2' => '5.7183527037985e-06',
        'eta'  => 0
    } );

for my $matrix_id ( 0..$#{ $matrix_product } ) {
    for my $row ( @{ $matrix_product->[$matrix_id] } ) {
        print( join( " ", @{ $row } ), "\n" );
    }
    print( "\n" ) if $matrix_id != $#{ $matrix_product };
}
