#!/bin/bash

export PERL5LIB=$(dirname "$0")/../../lib
export PDBX_FILE_1=$(dirname "$0")/../inputs/amino-acids/serine-001.cif
export PDBX_FILE_2=$(dirname "$0")/../inputs/amino-acids/aspartic-acid-001.cif

perl <<'END'
#!/usr/bin/perl

use strict;
use warnings;

use Data::Dumper;

use AtomSite;

$Data::Dumper::Sortkeys = 1;
$Data::Dumper::Indent = 1;

my $atom_site_1 = AtomSite->new();
my $atom_site_2 = AtomSite->new();

$atom_site_1->open( $ENV{PDBX_FILE_1} );
$atom_site_2->open( $ENV{PDBX_FILE_2} );

$atom_site_1->append( [ $atom_site_2 ] );

print Dumper $atom_site_1->{'_atoms'};

END
