package ForceField::Generate;

use strict;
use warnings;

use Exporter qw( import );
our @EXPORT_OK = qw( obtain_force_field_data
                     force_field_parameters );

use List::Util qw( uniq );

use PDBxParser qw( pdbx_loop_to_array
                   obtain_pdbx_line_new
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

sub force_field_parameters
{
    my ( $pdbx_file ) = @_;

    my $force_field_data = obtain_force_field_data( $pdbx_file );

    my %force_field_parameters = ();

    # Force field constants.
    # my $SOFT_EPSILON = $force_field_data->{'_[local]_force_field'}{'soft_epsilon'};
    # my $SOFT_N = $force_field_data->{'_force_field'}{'soft_n'};
    # my $LJ_K = $force_field_data->{'_force_field'}{'lj_k'};
    # my $C_K = $force_field_data->{'_force_field'}{'c_k'};
    # my $H_K = $force_field_data->{'_force_field'}{'h_k'};
    # my $T_K = $force_field_data->{'_force_field'}{'t_k'};
    # my $CUTOFF_ATOM = $force_field_data->{'_force_field'}{'cutoff_atom'};
    # my $CUTOFF_START = $force_field_data->{'_force_field'}{'cutoff_start'};
    # my $CUTOFF_END = $force_field_data->{'_force_field'}{'cutoff_end'};

    # use Data::Dumper; print STDERR Dumper $SOFT_EPSILON;

    return;
}

1;
