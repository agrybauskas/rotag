#!/bin/bash

export PERL5LIB=$(dirname "$0")/../../lib

perl <<'END'
#!/usr/bin/perl

use strict;
use warnings;

use Data::Dumper;

use Atom;

$Data::Dumper::Sortkeys = 1;
$Data::Dumper::Indent = 1;

my $atom = Atom->new( { 'id' => 1,
                        'type_symbol' => 'C',
                        'label_atom_id' => 'CA',
                        'Cartn_x' => 0.000,
                        'Cartn_y' => 0.000,
                        'Cartn_z' => 0.000 } );

print Dumper $atom;

END
