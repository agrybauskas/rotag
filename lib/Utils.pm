package Utils;

use strict;
use warnings;

use Exporter qw( import );
our @EXPORT_OK = qw( populate_multi_hash
                     retrieve_multi_hash );

sub populate_multi_hash
{
    my ( $hash, $keys ) = @_;
    for my $key ( @{ $keys } ) {

    }
}

sub retrieve_multi_hash
{
    my ( $hash, $keys, $options ) = @_;
    my ( $default ) = ( $options->{'default'} );
    if( $#{ $keys } ) {
        return defined $hash->{$keys->[0]} ?
               $hash->{$keys->[0]} :
               $options->{'default'};
    }
    return retrieve_multi_hash( $hash,
                                [ map { $keys->[$_] } 1..$#{ $keys } ],
                                $options );
}

1;
