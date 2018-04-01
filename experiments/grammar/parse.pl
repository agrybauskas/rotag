#!/usr/bin/perl

use strict;
use warnings;

use CommandLineParser;

my $input = "resname ASN";

my $parser = new CommandLineParser();
$parser->parser( $input );
