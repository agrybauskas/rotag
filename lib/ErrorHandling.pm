package ErrorHandling;

use strict;
use warnings;

use Exporter qw( import );
our @EXPORT_OK = qw( parse_errors
                     parse_warnings );

our $VERSION = '1.0.0';

# ------------------------------ Error handling ------------------------------- #

sub parse_errors
{
    my ( $args ) = @_;

    my ( $program, $filename, $message ) = (
        $args->{'program'},
        $args->{'filename'},
        $args->{'message'},
    );

    $program =~ s/^.+\/(\w+)$/$1/gsxm;

    print {*STDERR} "ERROR: $program: '$filename' - $message";

    exit 1;
}

sub parse_warnings
{
    my ( $args ) = @_;

    my ( $program, $filename, $message ) = (
        $args->{'program'},
        $args->{'filename'},
        $args->{'message'},
    );

    $program =~ s/^.+\/(\w+)$/$1/gsxm;

    print {*STDERR} "WARNING: $program: '$filename' - $message";

    return;
}

1;
