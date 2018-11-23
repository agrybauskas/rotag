package AtomSite;

use strict;
use warnings;

use Exporter qw( import );
our @EXPORT_OK = qw( append_atom_site
                     prepare_for_calc );

use BondProperties qw( hybridization );
use ConnectAtoms qw( connect_atoms );
use SidechainModels qw( rotation_only );

use Version qw( $VERSION );

our $VERSION = $VERSION;

# ------------------------ Manipulations with atom site ----------------------- #

#
# Appends atom site and checks if there are no duplicate ids.
# Input:
#     $atom_sites - array of atom site data structures (see PDBxParser.pm);
#     $options->{'renumber_from'} - renumbers atoms starting from certain int.
# Output:
#     %atom_site - concatinated atom site.
#

sub append_atom_site
{
    my ( $atom_sites, $options ) = @_;
    my ( $renumber_from ) = $options->{'renumber_from'};

    $renumber_from //= 0;

    my %atom_site;
    for my $atom_site ( @{ $atom_sites } ) {
        for my $atom_id ( sort keys %{ $atom_site } ) {
            if( ! exists $atom_site{$atom_id} ) {
                if( $renumber_from ) {
                    $atom_site{$renumber_from} = $atom_site->{$atom_id};
                    $atom_site{$renumber_from}{'id'} = $renumber_from;
                    $renumber_from++;
                } else {
                    $atom_site{$atom_id} = $atom_site->{$atom_id};
                }
            } else {
                die "atom site data structure with identical ids are being " .
                    "concatinated";
            }
        }
    }

    return \%atom_site;
}

#
# Prepares atom site data structure for calculations - adds connections,
# assigns hybridizations and rotational bonds.
# Input:
#     $atom_site - atom site data structure (see PDBxParser.pm).
# Output:
#     hybridization, connections and conformation attributes are added to the
#     $atom_site.
#

sub prepare_for_calc
{
    my ( $atom_site ) = @_;

    connect_atoms( $atom_site );
    hybridization( $atom_site );
    rotation_only( $atom_site );

    return;
}

1;
