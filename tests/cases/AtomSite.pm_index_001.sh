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

my $index_table = AtomSite->index( $atom_site );

for my $attribute ( sort keys %{ $index_table } ) {
    for my $value ( sort keys %{ $index_table->{$attribute} } ) {
        printf "%s %s %s\n",
               $attribute,
               $value,
               join( ',', sort @{ $index_table->{$attribute}{$value} } );
    }
}

END
