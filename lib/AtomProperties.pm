package AtomProperties;

use strict;
use warnings;

use Exporter qw( import );
our @EXPORT_OK = qw( sort_atom_names );

use Carp qw( confess );

use Version qw( $VERSION );

our $VERSION = $VERSION;

# -------------------------------- Atom names --------------------------------- #

#
# Sorts atom names by their hierarchical rules.
# Input:
#     $atom_names - list of atom names;
#     $options{sort_type} - sorts by atom type, greek letter, number or their
#     combinations.
# Output:
#     @sorted_names - sorted list of atom names.
#

sub sort_atom_names
{
    my ( $atom_names, $options ) = @_;
    my ( $sort_type ) = ( $options->{'sort_type'} );

    $sort_type //= 'tgn';

    # First priority is decided by atom type: S > P > O > N > C > H.
    # Second priority - by greek letter: A > B > G > D > E > Z > H.
    # TODO: look for more greek letters that are in PDBx.
    # TODO: hetatoms.
    # TODO: other atoms.
    # Third priority - by numeration: 1 > 2 > 3 and etc.
    # This priority is achieved by assinging first, second and third priorities
    # to numbers. Then iteratively is sorted by priorities.
    my %atom_type_priority =
        ( 'H' => 1, 'X' => 2, 'C' => 3, 'N' => 4, 'O' => 5, 'P' => 6, 'S' => 7,
          'MG' => 8 );
    my %greek_letter_priority =
        ( 'H' => 1, 'Z' => 2, 'E' => 3, 'D' => 4, 'G' => 5, 'B' => 6, 'A' => 7,
          q() => 8, );

    my @atom_types =
        sort { $atom_type_priority{$a} <=> $atom_type_priority{$b} }
        keys %atom_type_priority;
    my @greek_letters =
        sort { $greek_letter_priority{$a} <=> $greek_letter_priority{$b} }
        keys %greek_letter_priority;

    # Decomposes each atom name by its components.
    my %atom_names;
    for my $atom_name ( @{ $atom_names } ) {
        my ( $atom_type ) =
            $atom_name =~ /(^\p{IsAlphabetic})\p{IsAlphabetic}?\d?/smx;
        my ( $greek_letter ) =
            $atom_name =~ /^\p{IsAlphabetic}(\p{IsAlphabetic}?)\d?/smx;
        my ( $number ) =
            $atom_name =~ /^\p{IsAlphabetic}\p{IsAlphabetic}?(\d?)/smx;
        $atom_names{$atom_name}{'type'} =
            $atom_type_priority{$atom_type};
        $atom_names{$atom_name}{'greek_letter'} =
            $greek_letter_priority{$greek_letter};
        $atom_names{$atom_name}{'number'} = $number;
    }

    # Sorts by rules of described in %atom_names.
    my @sorted_names;
    if( $sort_type eq 'tgn') { # By type, then greek letter and then number.
        @sorted_names =
          sort {
            $atom_names{$b}{'type'} <=> $atom_names{$a}{'type'} ||
            $atom_names{$b}{'greek_letter'} <=> $atom_names{$a}{'greek_letter'}||
            $atom_names{$a}{'number'} cmp $atom_names{$b}{'number'} }
              @{ $atom_names };
    } elsif( $sort_type eq 'gn' ) { # By greek letter and then number.
        @sorted_names =
          sort {
            $atom_names{$b}{'greek_letter'} <=> $atom_names{$a}{'greek_letter'}||
            $atom_names{$a}{'number'} cmp $atom_names{$b}{'number'} }
              @{ $atom_names };
    } else {
        confess "there is no such '${sort_type}' option when sorting atom names";
    }

    return \@sorted_names;
}

1;
