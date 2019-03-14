#!/usr/bin/perl

use strict;
use warnings;

use Data::Dumper;

use ForceField::Parameters;

$Data::Dumper::Sortkeys = 1;
$Data::Dumper::Indent = 1;

my $PARAMETERS = Parameters->new();

print Dumper $PARAMETERS;
