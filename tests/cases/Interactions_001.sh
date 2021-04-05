#!/bin/bash

export PDBx_file=$(dirname "$0")/../inputs/amino-acids/aspartic-acid-001.cif

perl <<'END'
use strict;
use warnings;

use Data::Dumper;

use PDBxParser qw( raw2indexed
                   obtain_pdbx_data );
use Interactions;

$Data::Dumper::Sortkeys = 1;
$Data::Dumper::Indent = 1;

my $pdbx_data = obtain_pdbx_data( "$ENV{PDBx_file}", [ '_atom_site' ] );
raw2indexed( $pdbx_data,
             { 'attributes' => { '_atom_site' => [ 'id' ] } } );
my $atom_site = $pdbx_data->{'_atom_site'}{'data'};

my $interactions = Interactions->new;
$interactions->add( $pdbx_data );
END
