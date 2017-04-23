#!/usr/bin/perl
  
use strict;
use warnings;

#
# Trying cos(x)^2 + sin(x)^2 = 1 type expression.
#

foreach ( ( 1..10 ) ) {
    system( "ginsh", qq( cos(x) + cos(x); ) );
}
