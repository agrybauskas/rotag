#!/bin/bash

export PERL5LIB=$(dirname "$0")/../../lib
export PDBX_FILE=$(dirname "$0")/../inputs/5svd.cif

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

my $filtered_atom_site =
   $atom_site->filter( { 'include' => { 'label_atom_id' => [ 'CA', 'CB' ] } } );

print Dumper $filtered_atom_site->extract( { 'data' => [ 'label_atom_id' ] } );

END
