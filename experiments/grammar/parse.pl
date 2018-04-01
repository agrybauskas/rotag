#!/usr/bin/perl

use strict;
use warnings;

use Data::Dumper;

use CommandLineParser;

my $input = "resid 10";

my $parser = new CommandLineParser();
$parser->parser( $input );

print Dumper $parser->{'include'}{'resid'};
