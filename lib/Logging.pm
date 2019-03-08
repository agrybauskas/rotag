package Logging;

use strict;
use warnings;

use Exporter qw( import );
our @EXPORT_OK = qw( error
                     info
                     warning );

use Version qw( $VERSION );

our $VERSION = $VERSION;

# ------------------------------ Error handling ------------------------------- #

sub error
{
    my ( $args ) = @_;

    my ( $program, $filename, $type, $message ) = (
        $args->{'program'},
        $args->{'filename'},
        $args->{'type'},
        $args->{'message'},
    );

    $type //= 'Error';

    if( defined $filename && $filename eq '-' ) { $filename = 'STDIN'; }

    $program =~ s/^.+\/(\w+)$/$1/gsxm;

    if( defined $filename ) {
        print {*STDERR} "$type: $program: $filename - $message";
    } else {
        print {*STDERR} "$type: $program: $message";
    }

    exit 1;
}

sub info
{
    my ( $args ) = @_;

    my ( $type, $message ) = (
        $args->{'type'},
        $args->{'message'},
    );

    $type //= 'Info';

    print {*STDERR} "$type: $message";

    return;
}

sub warning
{
    my ( $args ) = @_;

    my ( $program, $filename, $type, $message ) = (
        $args->{'program'},
        $args->{'filename'},
        $args->{'type'},
        $args->{'message'},
    );

    $type //= 'Warning';

    if( defined $filename && $filename eq '-' ) { $filename = 'STDIN'; }

    if( defined $filename ) {
        print {*STDERR} "$type: $program: $filename - $message";
    } elsif( defined $program ) {
        $program =~ s/^.+\/(\w+)$/$1/gsxm;
        print {*STDERR} "$type: $program: $message";
    } else {
        print {*STDERR} "$type: $message";
    }

    return;
}

1;
