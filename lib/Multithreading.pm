package Multithreading;

use strict;
use warnings;

use Exporter qw( import );
our @EXPORT_OK = qw( multithreading
                     forking );

use Carp;
use Parallel::ForkManager;
use threads;

use Version qw( $VERSION );

our $VERSION = $VERSION;

# -------------------- Preparing data for multi-threading --------------------- #

#
# Divides data arrays into smaller arrays that are being passed to multithreading
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
        my @array = @{ $array };

        # Arrays have to have equal lengths.
        if( defined $array_length && scalar( @array ) ne $array_length ) {
            confess 'list of arrays have different lengths.';
        }
        $array_length =  scalar @array;

        # Splits the array into blocks/chunks.
        my @array_blocks;
        my $max_block_size = int( scalar @array / $threads );

        # If block size is smaller that the number of threads, then
        # thread number is reduced.
        my $reduce_threads = 0;
        if( ! $max_block_size ) {
            $reduce_threads = scalar @array;
            $max_block_size = 1;
        }

        for my $i ( 0..$threads-$reduce_threads-1 ) {
            my $block_start = $i * $max_block_size;
            my $block_end = $block_start + $max_block_size - 1;
            if( $i ne $threads-$reduce_threads-1 ) {
                push @array_blocks, [ @array[$block_start..$block_end] ];
            } else {
                $block_end = $#array;
                push @array_blocks, [ @array[$block_start..$block_end] ];
            }
        }

        push @list_of_array_blocks, \@array_blocks;
    }

    return \@list_of_array_blocks;
}

# ----------------------------- Multi-threading ------------------------------- #

#
# Performs multithreading where function, arguments and data blocks are given.
# Input:
#     $function - reference to the function that will be run;
#     $arguments - arguments that will be passed to a function;
#     $divisible_arrays - array of data;
#     $threads - number of threads.
# Output:
#     @joined_block_results - results from all threads.
#

sub multithreading
{
    my ( $function, $arguments, $divisible_arrays, $threads ) = @_;

    my $allowed_array_blocks =
        divide_arrays_into_blocks( $divisible_arrays, $threads );

    my @block_results;
    for my $i ( 0..$threads-1 ) {
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
# Performs forking (similar idea to multithreading()) where function,
# arguments and data blocks are given.
# Input:
#     $function - reference to the function that will be run;
#     $arguments - arguments that will be passed to a function;
#     $divisible_arrays - array of data;
#     $proc_num - number of processes that will be run simultaneously.
# Output:
#     @joined_block_results - results from all threads.
#

sub forking
{
    my ( $function, $arguments, $divisible_arrays, $proc_num )= @_;
}

1;
