package Utils;

use strict;
use warnings;

use lib "../../lib/perl";
use CifParser;

sub select_atom_data
{
    my $attribute_selector = shift;
    my $data_selector = shift;
    my @cif = @_;

    # Parses selector argument from string  to proper array.
    my %attribute_selector = ( map { $_->[0] => [ split( ",", @$_[1] ) ] }
                               map { [ split( " ", $_ ) ] }
                               split( "&", $attribute_selector ) );
    my @data_selector = split( ",", $data_selector );

    # Selects atoms for further analysis.
    my @selected_atom_data =
	@{ CifParser::select_atom_data(
	       \@data_selector,
	       &CifParser::filter_atoms(
		   \%attribute_selector,
		   &CifParser::obtain_atom_site( @_ ) ) ) };

    return \@selected_atom_data;
}

1;
