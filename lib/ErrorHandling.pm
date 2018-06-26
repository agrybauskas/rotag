package ErrorHandling;

use strict;
use warnings;

use Exporter qw( import );
our @EXPORT_OK = qw( parse_errors
                     parse_warnings );

# ------------------------------ Error handling ------------------------------- #

sub parse_errors
{
    my ( $args ) = @_;

    my ( $program, $filename, $message ) = (
        $args->{'program'},
        $args->{'filename'},
        $args->{'message'}
    );

    $program =~ s/^.+\/(\w+)$/$1/g;

    print STDERR "ERROR: $program: '$filename' - $message";
    exit 1;
}

sub parse_warnings
{
    my ( $args ) = @_;

    my ( $program, $filename, $message ) = (
        $args->{'program'},
        $args->{'filename'},
        $args->{'message'}
    );

    $program =~ s/^.+\/(\w+)$/$1/g;

    print STDERR "WARNING: $program: '$filename' - $message";
}
