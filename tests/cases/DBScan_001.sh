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

$cluster->dbscan( 1.0 );

for my $id ( sort { $a <=> $b } @{ $cluster->get_ids() } ) {
    print $id . " " . join( ',', @{ $cluster->get_point( $id )->{'coord'} } ) .
          " " . $cluster->get_point( $id )->{'label'} . "\n";
}
