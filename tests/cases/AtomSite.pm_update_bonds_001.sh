#!/bin/bash

export PERL5LIB=$(dirname "$0")/../../lib
export PDBX_FILE=$(dirname "$0")/../inputs/amino-acids/serine-001.cif

perl <<'END'
#!/usr/bin/perl

use strict;
use warnings;

use Data::Dumper;

use AtomSite;

$Data::Dumper::Sortkeys = 1;
$Data::Dumper::Indent = 1;

my $atom_site = AtomSite->new();
$atom_site->open( $ENV{PDBX_FILE} );
$atom_site->update_bonds();

print Dumper $atom_site->{'_atoms'}

END
