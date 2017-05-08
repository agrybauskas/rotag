package Utils;

use strict;
use warnings;

use Exporter qw( import );
our @EXPORT_OK = qw( select_atom_data );

use lib "../../lib/perl";
use CifParser qw( filter_atoms
                  obtain_atom_site );

sub select_atom_data
{
    my ( $attribute_selector, $data_selector, $cif_file ) = @_;

    # Parses selector argument from string  to proper array.
    my %attribute_selector = ( map { $_->[0] => [ split( ",", @$_[1] ) ] }
                               map { [ split( " ", $_ ) ] }
                               split( "&", $attribute_selector ) );

    my @data_selector = split( ",", $data_selector );

    # Selects atoms for further analysis.
    my @selected_atom_data =
	@{ CifParser::select_atom_data(
	       \@data_selector,
	       &filter_atoms( \%attribute_selector,
			      &obtain_atom_site( $cif_file ) ) ) };

    return \@selected_atom_data;
}

1;
