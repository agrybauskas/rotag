package Multiprocessing;

use strict;
use warnings;

use Exporter qw( import );
our @EXPORT_OK = qw( divide_arrays_into_blocks
                     forking
                     threading );

use Carp;
use Parallel::ForkManager;
use POSIX qw( ceil );
use threads;

use Version qw( $VERSION );

our $VERSION = $VERSION;

# ---------------------- Preparing data for threading ------------------------- #

#
# Divides data arrays into smaller arrays that are being passed to threading
# function.
# Input:
#     $arrays - list of arrays of data;
#     $threads - number of threads.
# Output:
#     @list_of_array_blocks - list of arrays of arrays of data.
#

sub divide_arrays_into_blocks
{
    my ( $arrays, $threads ) = @_;

    my $array_length;
    my @list_of_array_blocks;
    for my $array ( @{ $arrays } ) {
        my ( @array ) = @{ $array };

        # Arrays have to have equal lengths.
        if( defined $array_length && scalar( @array ) ne $array_length ) {
            confess 'list of arrays have different lengths.';
        }

        $array_length =  scalar @array;

        # Splits the array into blocks/chunks.
        my $max_block_size = ceil( $array_length / $threads );

        my @array_blocks;
        for( my $i = 0; $i < $array_length; $i += $max_block_size ) {
            my $block_start = $i;
            my $block_end = $i + $max_block_size - 1;

            if( $block_end >= $array_length ) {
                $block_end = $array_length - 1;
            }

            push @array_blocks, [ @array[$block_start..$block_end] ];
        }

        push @list_of_array_blocks, \@array_blocks;
    }

    return \@list_of_array_blocks;
}

# -------------------------------- Threading ---------------------------------- #

#
# Performs threading where function, arguments and data blocks are given.
# Input:
#     $function - reference to the function that will be run;
#     $arguments - arguments that will be passed to a function;
#     $divisible_arrays - array of data;
#     $threads - number of threads.
# Output:
#     @joined_block_results - results from all threads.
#

sub threading
{
    my ( $function, $arguments, $divisible_arrays, $threads ) = @_;

    my $allowed_array_blocks =
        divide_arrays_into_blocks( $divisible_arrays, $threads );

    my @block_results;
    for my $i ( 0..$#{ $allowed_array_blocks->[0] } ) {
        my $thread_task =
            threads->create( $function,
                             $arguments,
                             [ map { $allowed_array_blocks->[$_][$i] }
                                   0..$#{$allowed_array_blocks} ] );
        push @block_results, $thread_task;
    }

    my @joined_block_results;
    for my $block_result ( @block_results ) {
        $block_result = $block_result->join();
        for( my $i = 0; $i <= $#{ $block_result }; $i++ ) {
            push @{ $joined_block_results[$i] }, @{ $block_result->[$i] };
        }
    }

    return \@joined_block_results;
}

# --------------------------------- Forking ----------------------------------- #

#
# Performs forking (similar idea to threading()) where function, arguments and
# data blocks are given.
# Input:
#     $function - reference to the function that will be run;
#     $arguments - arguments that will be passed to a function;
#     $data - array of single/multiple data;
#     $threads - number of processes that will be run simultaneously.
# Output:
#     @results - results from all processes.
#

sub forking
{
    my ( $function, $arguments, $data, $threads )= @_;

    my $pm = Parallel::ForkManager->new( $threads );

    my @results;

    # Rule to retrieve data.
    $pm->run_on_finish( sub {
         my ( $pid, $exit_code, $ident, $exit_signal,
              $core_dump, $data_structure_reference ) = @_;
         for my $i ( 0..$#{ $$data_structure_reference } ) {
             push @{ $results[$i] }, @{ $$data_structure_reference->[$i] };
         }
    } );

  DATA:
    for my $i ( 0..$#{ $data->[0] } ) {
        $pm->start and next DATA;
        my $array = $function->( $arguments, $data->[0][$i], $data->[1][$i] );
        $pm->finish( 0, \$array );
    }
    $pm->wait_all_children;

    return \@results;
}

1;
