package ForceField::Generate;

use strict;
use warnings;

use Exporter qw( import );
our @EXPORT_OK = qw( obtain_force_field_data );

use List::Util qw( uniq );

use PDBxParser qw( obtain_pdbx_line_new
                   obtain_pdbx_loop );

use Version qw( $VERSION );

our $VERSION = $VERSION;

sub obtain_force_field_data
{
    my ( $pdbx_file ) = @_;

    my @looped_categories;
    my @unlooped_items;

    my $is_looped = 0;

    local @ARGV = ( $pdbx_file );
    while( <> ) {
        if( m/^loop_/ ) {
            $is_looped = 1;
        } elsif( $is_looped && m/^(_\S+)\./ ) {
            my $category = $1;
            push @looped_categories, $category;
        } elsif( m/(^_\S+\.\S+)/ ) {
            my $item = $1;
            push @unlooped_items, $item;
        } else {
            $is_looped = 0;
        }
    }

    @looped_categories = uniq @looped_categories;

    my $pdbx_line_data = obtain_pdbx_line_new( $pdbx_file, \@unlooped_items );
    my $pdbx_loop_data = obtain_pdbx_loop( $pdbx_file, \@looped_categories );

    return { %{ $pdbx_line_data }, %{ $pdbx_loop_data } };
}

1;
