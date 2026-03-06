#!/usr/bin/perl

use strict;
use warnings;

use Data::Dumper;

use DBScan;

$Data::Dumper::Sortkeys = 1;
$Data::Dumper::Indent = 1;

my $cluster = DBScan->new( [ 1.0, 2.0, 3.0, 5.0, 6.0 ] );

print Dumper $cluster;
