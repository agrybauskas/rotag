#!/usr/bin/perl

use strict;
use warnings;

use Data::Dumper;

use CommandLineParser;

my $input = 'resname ASP,ASN';

my $parser = new CommandLineParser();
$parser->parser( $input );

print Dumper $parser->{'include'}{'resname'};
