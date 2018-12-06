#!/bin/bash

export PERL5LIB=$(dirname "$0")/../../lib

perl <<'END' 2>&1 | sed 's/0x[0-9a-f]\{12\}/<hex>/g'
#!/usr/bin/perl

use strict;
use warnings;

use Data::Dumper;

use Atom;

$Data::Dumper::Sortkeys = 1;
$Data::Dumper::Indent = 1;

my $atom = Atom->new( [ 'id' => 1 ] );

print Dumper $atom;

END
