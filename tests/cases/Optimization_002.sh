#!/usr/bin/perl

use strict;
use warnings;

use Data::Dumper;

use ForceField::Optimization;
use ForceField::Parameters;

$Data::Dumper::Sortkeys = 1;
$Data::Dumper::Indent = 1;

my $parameters = Parameters->new();

my $optimization = Optimization->new();
$optimization->parameters( $parameters );

print Dumper $optimization;
