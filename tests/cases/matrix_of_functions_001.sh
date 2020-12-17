#!/usr/bin/perl

use strict;
use warnings;

use Data::Dumper;

use LinearAlgebra qw( matrix_of_functions );

$Data::Dumper::Sortkeys = 1;
$Data::Dumper::Indent = 1;

print Dumper matrix_of_functions( \&cos, 3, 3 );
