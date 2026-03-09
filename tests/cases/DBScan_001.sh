#!/usr/bin/perl

use strict;
use warnings;

use Data::Dumper;

use DBScan;

$Data::Dumper::Sortkeys = 1;
$Data::Dumper::Indent = 1;

my $cluster =
    DBScan->new( [ [ 1.0 ], [ 2.0 ], [ 3.0 ], [ 5.0 ], [ 10.0 ], [ 11.0 ],
                   [ 20.0 ] ] );

print Dumper $cluster;

$cluster->dbscan( 1.0 );

print Dumper $cluster;
