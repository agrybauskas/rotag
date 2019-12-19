#!/usr/bin/perl

use strict;
use warnings;

use Data::Dumper;

use ForceField::Optimization;

$Data::Dumper::Sortkeys = 1;
$Data::Dumper::Indent = 1;

my $optimization = Optimization->new();

print Dumper $optimization;
