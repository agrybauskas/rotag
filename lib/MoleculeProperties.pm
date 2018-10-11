package MoleculeProperties;

use strict;
use warnings;

use Exporter qw( import );
our @EXPORT_OK = qw( %RESIDUE_ATOMS );

use Version qw( $VERSION );

our $VERSION = $VERSION;

my %RESIDUE_ATOMS = (
    '' => {
        'mandatory' =>
            [

            ],
        'optional'  =>
            [

            ],
    },
);

1;
