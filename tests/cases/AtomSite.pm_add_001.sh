#!/bin/bash

export PERL5LIB=$(dirname "$0")/../../lib
export PDBX_FILE_1=$(dirname "$0")/../inputs/amino-acids/serine-001.cif
export PDBX_FILE_2=$(dirname "$0")/../inputs/amino-acids/aspartic-acid-001.cif

perl <<'END'
#!/usr/bin/perl

use strict;
use warnings;

use Data::Dumper;

use Atom;
use AtomSite;

$Data::Dumper::Sortkeys = 1;
$Data::Dumper::Indent = 1;

my $atom_site = AtomSite->new();

$atom_site->add(
    [ Atom->new( { 'id' => 3,
                 'type_symbol' => 'C',
                 'label_atom_id' => 'CA',
                 'Cartn_x' => 0.000,
                 'Cartn_y' => 0.000,
                 'Cartn_z' => 0.000 } ) ] );

print Dumper $atom_site;

END
