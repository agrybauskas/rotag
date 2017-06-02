package PseudoAtoms;

use strict;
use warnings;

use Exporter qw( import );
our @EXPORT_OK = qw( generate_pseudo );

use lib qw( ./ );
use Data::Dumper;
my $parameter_file = "../../parameters/vdw_radii.csv";

# --------------------------- Generation of pseudo-atoms ---------------------- #

#
# Generates pseudo-atoms from side chain models that are written in equation
# form.
#

sub generate_pseudo
{
    my ( $atom_site ) = @_;

    print Dumper $atom_site;
}
