#!/bin/bash

export PERL5LIB=$(dirname "$0")/../../lib
export PDBX_FILE_1=$(dirname "$0")/../inputs/5svd.cif
export PDBX_FILE_2=$(dirname "$0")/../inputs/1bvi.cif

perl <<'END'
#!/usr/bin/perl

use strict;
use warnings;

use Data::Dumper;

use AtomSite;

$Data::Dumper::Sortkeys = 1;
$Data::Dumper::Indent = 1;

my $atom_site = AtomSite->new();
my $atom_site_1 = AtomSite->new();
my $atom_site_2 = AtomSite->new();

$atom_site_1->open( $ENV{PDBX_FILE_1} );
$atom_site_2->open( $ENV{PDBX_FILE_2} );

$atom_site->append( [ $atom_site_1->{'atoms'},
                      $atom_site_2->{'atoms'} ], 'r' );

print Dumper $atom_site;

END
