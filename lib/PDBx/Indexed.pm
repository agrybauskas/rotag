package PDBx::Indexed;

use strict;
use warnings;

use Carp;

use PDBx::Raw;

# ------------------------- Constructors/Destructors -------------------------- #

sub new
{
    my ( $class, $pdbx_file, $data_identifier, $options ) = @_;
    my $pdbx_raw = PDBx::Raw->new( $pdbx_file, $data_identifier );
    my $self = indexed( $pdbx_raw, $options );
    return bless $self, $class;
}

# ----------------------------- Setters/Getters ------------------------------- #

sub get_category
{
    my ( $self, $category ) = @_;
    if( exists $self->{$category} ) {
        return $self->{$category};
    } else {
        return {};
    }
}

sub set_category
{
    my ( $self, $category, $data ) = @_;
    $self->{$category} = $data;
}

# --------------------------------- Methods ----------------------------------- #

sub read
{
    my ( $self, $pdbx_file, $data_identifier, $options ) = @_;
    my $pdbx_raw = PDBx::Raw->new( $pdbx_file, $data_identifier );
    my $pdbx_data_indexed = indexed( $pdbx_raw, $options );
    for my $category ( keys %{ $pdbx_data_indexed } ) {
        $self->{$category} = $pdbx_data_indexed->{$category};
    }
}

#
# Takes PDBx and indexes them by defined attributes.
# Input:
#     $attributes - combination of attribute data that serves as unique key.
#     $options->{'read_until_end'} - reads whole pdbx file or stdin.
# Output:
#     %indexed - returns indexed pdbx;
#

sub indexed
{
    my ( $self, $options ) = @_;
    my ( $attributes, $read_until_end ) = (
        $options->{'attributes'},
        $options->{'read_until_end'}
    );

    $attributes //= {} unless $attributes;
    $read_until_end //= 0;

    my %attributes = %{ $attributes };
    my $pdbx_raw = $self;

    if( ! defined $pdbx_raw || ! %{ $pdbx_raw } ) {
        return {};
    }

    my %indexed = ();

    for my $category ( keys %{ $pdbx_raw } ) {
        my $keys;
        if( ! exists $attributes->{$category} ) {
            $keys = [ 'auto_increment' ]; # Special key.
        } else {
            $keys = $attributes->{$category};
        }

        my @attributes = @{ $pdbx_raw->{$category}{'metadata'}{'attributes'} };
        my $is_unique = $pdbx_raw->{$category}{'metadata'}{'is_unique'};
        $is_unique //= 1;
        my @data = @{ $pdbx_raw->{$category}{'data'} };

        # Creates special data structure.
        my @data_row;
        my %data_row;

        my $attribute_count = scalar @attributes;
        my $data_count = scalar @data;

        # Determines the positions of unique keys in attribute list.
        my @attribute_pos;
        for my $key ( @{ $keys } ) {
            for( my $i = 0; $i <= $#attributes; $i++ ) {
                if( $attributes[$i] eq $key ) {
                    push @attribute_pos, $i;
                    last;
                }
            }
        }

        my $auto_increment = 0;
        my %current_indexed = ();
        $current_indexed{'metadata'}{'is_loop'} =
            $pdbx_raw->{$category}{'metadata'}{'is_loop'};
        for( my $pos = 0; $pos < $data_count - 1; $pos += $attribute_count ) {
            @data_row = @{ data[$pos..$pos+$attribute_count-1] };
            %data_row = ();
            for( my $col = 0; $col <= $#data_row; $col++ ) {
                $data_row{$attributes[$col]} = $data_row[$col];
            }

            my $current_key = join q{,}, map { $data_row[$_] } @attribute_pos;
            $current_key = $auto_increment unless $current_key;

            if( ! exists $current_indexed{'data'}{$current_key} ) {
                if( ! $is_unique ) {
                    $current_indexed{'data'}{$current_key} = [ { %data_row } ];
                } else {
                    $current_indexed{'data'}{$current_key} = { %data_row };
                }
            } else {
                if( ! $is_unique ) {
                    push @{ $current_indexed{'data'}{$current_key} },
                          { %data_row };
                } else {
                    confess 'unique key supposed to be unique.';
                }
            }

            $auto_increment++;
        }

        $indexed{$category} = { %current_indexed };
    }

    return \%indexed;
}

sub to_pdbx
{
    my ( $self ) = @_;
    bless( PDBx::Raw::raw( $self ), 'PDBx::Raw' )->to_pdbx();
}

sub to_csv
{
    my ( $self, $options ) = @_;
    bless( PDBx::Raw::raw( $self ), 'PDBx::Raw' )->to_csv( $options );
}

1;
