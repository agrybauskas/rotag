package Utils;

use strict;
use warnings;

use lib "../../lib/perl";
use CifParser;

sub select_atom_data
{
    my $attribute_filter = shift;
    my $attribute_select = shift;
    my @cif = @_;

    # Parses selector argument from string  to proper array.
    my @attribute_filter = [ map { $_->[0] => [ split( ",", @$_[1] ) ] }
			     map { [ split( " ", $_ ) ] }
			     split( "&", $attribute_filter ) ];

    my @attribute_select = [ split( ",", $attribute_select ) ];

    # Selects atoms for further analysis.
    my @selected_atom_data =
	@{ CifParser::select_atom_data( @attribute_filter,
					@attribute_select,
					@cif ) };

    return \@selected_atom_data;
}

1;
