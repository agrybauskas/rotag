#!/usr/bin/perl

use strict;
use warnings;

use Data::Dumper;

use DBScan;

$Data::Dumper::Sortkeys = 1;
$Data::Dumper::Indent = 1;

my $cluster =
    DBScan->new( [ [ -31.184, -6.586, -1.238 ],
                   [ -31.549400, -8.643799, 6.083401 ],
                   [ -30.882276, -8.018614, 5.698467 ],
                   [ -30.770998, -7.217951, 5.123454 ],
                   [ -31.258069, -6.547636, 4.577997 ],
                   [ -32.157446, -6.263706, 4.270443 ] ] );

$cluster->dbscan( 0.993 );

for my $id ( sort { $a <=> $b } @{ $cluster->get_ids() } ) {
    print $id . " " . join( ',', @{ $cluster->get_point( $id )->{'coord'} } ) .
          " " . $cluster->get_point( $id )->{'label'} . "\n";
}
