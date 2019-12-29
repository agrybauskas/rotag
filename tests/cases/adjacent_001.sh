#!/usr/bin/perl

use strict;
use warnings;

use Data::Dumper;

use Combinatorics qw( adjacent );

$Data::Dumper::Sortkeys = 1;
$Data::Dumper::Indent = 1;

print Dumper adjacent( [ [ 1, 2 ], [ 3, 4 ] ] );
