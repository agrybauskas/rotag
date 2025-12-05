package StructConn;

use strict;
use warnings;

use Exporter qw( import );
our @EXPORT_OK = qw( create_struct_conn );

use Carp;

use Version qw( $VERSION );

our $VERSION = $VERSION;

# ----------------------------------------------------------------------------- #

sub create_struct_conn
{
    my ( $atom1, $atom2, $type, $id ) = @_;
}

1;
