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
print Dumper create( { 'id' => 1,
                       'type_symbol' => 'C',
                       'label_atom_id' => 'CA',
                       'label_alt_id' => '.',
                       'label_comp_id' => 'SER',
                       'label_asym_id' => 'A',
                       'label_entity_id' => '1',
                       'label_seq_id' => '1',
                       'Cartn_x' => 0.0,
                       'Cartn_y' => 0.0,
                       'Cartn_z' => 0.0,
                       'pdbx_PDB_model_num' => 1 } );

END
