package ForceField::Generate;

use strict;
use warnings;

use Exporter qw( import );
our @EXPORT_OK = qw( obtain_force_field_data );

use List::Util qw( uniq );

use Version qw( $VERSION );

our $VERSION = $VERSION;

sub obtain_force_field_data
{
    my ( $pdbx_file ) = @_;

    my @lined_categories;
    my @looped_categories;

    my $is_looped = 0;

    local @ARGV = ( $pdbx_file );
    while( <> ) {
        if( m/^loop_/ ) {
            $is_looped = 1;
        } elsif( $is_looped && m/^_(\S+)\./ ) {
            my $category = $1;
            push @looped_categories, $category;
        } elsif( m/^_(\S+)\./ ) {
            my $category = $1;
            push @lined_categories, $category;
        } else {
            $is_looped = 0;
        }
    }

    @lined_categories = uniq @lined_categories;
    @looped_categories = uniq @looped_categories;
}

1;
