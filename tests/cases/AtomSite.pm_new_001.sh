#!/bin/bash

export PERL5LIB=$(dirname "$0")/../../lib

perl <<'END'
#!/usr/bin/perl

use strict;
use warnings;

use Data::Dumper;

use AtomSite;

$Data::Dumper::Sortkeys = 1;
$Data::Dumper::Indent = 1;

my $atom_site = AtomSite->new();

print Dumper $atom_site;

END
