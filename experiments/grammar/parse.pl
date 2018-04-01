#!/usr/bin/perl

use strict;
use warnings;

use Data::Dumper;

use CommandLineParser;

my $input = '! resid 20 | resname ASP,ASN | resid 2-10,1';

my $parser = new CommandLineParser();
$parser->parser( $input );

print Dumper $parser->{'include'}{'resname'};
print Dumper $parser->{'include'}{'resid'};
