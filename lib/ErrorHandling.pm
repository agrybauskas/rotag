package ErrorHandling;

use strict;
use warnings;

use Exporter qw( import );
our @EXPORT_OK = qw( parse_errors
                     parse_warnings );

use Version qw( $VERSION );

our $VERSION = $VERSION;

# ------------------------------ Error handling ------------------------------- #

sub parse_errors
{
    my ( $args ) = @_;

    my ( $program, $filename, $type, $message ) = (
        $args->{'program'},
        $args->{'filename'},
        $args->{'type'},
        $args->{'message'},
    );

    $type //= 'ERROR';

    if( $filename eq '-' ) { $filename = 'STDIN'; }

    $program =~ s/^.+\/(\w+)$/$1/gsxm;

    print {*STDERR} "$type: $program: $filename - $message";

    exit 1;
}

sub parse_warnings
{
    my ( $args ) = @_;

    my ( $program, $type, $filename, $message ) = (
        $args->{'program'},
        $args->{'filename'},
        $args->{'type'},
        $args->{'message'},
    );

    $type //= 'WARNING';

    if( $filename eq '-' ) { $filename = 'STDIN'; }

    $program =~ s/^.+\/(\w+)$/$1/gsxm;

    print {*STDERR} "$type: $program: '$filename' - $message";

    return;
}

1;
