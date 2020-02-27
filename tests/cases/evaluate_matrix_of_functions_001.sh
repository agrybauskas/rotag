#!/usr/bin/perl

use strict;
use warnings;

use Data::Dumper;

use ForceField::Parameters;
use LinearAlgebra qw( evaluate_matrix_of_functions );

$Data::Dumper::Sortkeys = 1;
$Data::Dumper::Indent = 1;

my $parameters = Parameters->new();
my $pi = $parameters->{'_[local]_constants'}{'pi'};

print Dumper evaluate_matrix_of_functions(
    [ [ \&CORE::cos, \&CORE::sin, \&CORE::cos ],
      [ \&CORE::sin, \&CORE::cos, \&CORE::sin ] ],
    [ [ 0, $pi, 0  ],
      [ $pi, 0, $pi ] ]
);
