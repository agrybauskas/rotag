package PDB;

use strict;
use warnings;

use Exporter qw( import );
our @EXPORT_OK = qw( obtain_pdb_atom_site );

use Version qw( $VERSION );

our $VERSION = $VERSION;

# --------------------------------- PDB parser -------------------------------- #

#
# Obtains pdb data and converts to atom site data structure.
# Input:
#     $pdb_file - PDB file path.
# Output:
#     %atom_site - atom site data structure.
#

sub obtain_pdb_atom_site
{
    my ( $pdbx_file ) = @_;

    my %atom_site = ();

    my $pdbx_PDB_model_num = 1;

    local @ARGV = ( $pdbx_file );
    while( <> ) {
        if( /^MODEL/ ) {
            my ( $model_id ) = m/^.{6}.{4}(.{4})$/;
            $pdbx_PDB_model_num = $model_id;
        } elsif( /^ATOM|^HETATM/ ) {
            # TODO: look if more '?!:' should be included in regex.
            my @atom_data = m/^(.{6})(.{5}).{1}(.{4})(.{1})(.{3}).{1}(.{1})(.{4}).{4}(.{8})(.{8})(.{8})(.{6})(.{6}).{10}(.{2})(?!:.{2})/;

            for my $atom_data_item ( @atom_data ) {
                $atom_data_item =~ s/ //g;
            }

            my ( $group_pdb, $id, $label_atom_id, $label_alt_id, $label_comp_id,
                 $label_asym_id, $label_seq_id, $cartn_x, $cartn_y, $cartn_z,
                 $occupancy, $B_iso_or_equiv_esd, $type_symbol,
                 $pdbx_formal_charge ) = @atom_data;

            $atom_site{'data'}{$id}{'group_PDB'} = $group_pdb;
            $atom_site{'data'}{$id}{'id'} = $id;
            $atom_site{'data'}{$id}{'type_symbol'} = $type_symbol;
            $atom_site{'data'}{$id}{'label_atom_id'} = $label_atom_id;
            $atom_site{'data'}{$id}{'label_alt_id'} = '.' unless $label_alt_id;
            $atom_site{'data'}{$id}{'label_comp_id'} = $label_comp_id;
            $atom_site{'data'}{$id}{'label_asym_id'} = $label_asym_id;
            $atom_site{'data'}{$id}{'label_seq_id'} = $label_seq_id;
            $atom_site{'data'}{$id}{'Cartn_x'} = $cartn_x;
            $atom_site{'data'}{$id}{'Cartn_y'} = $cartn_y;
            $atom_site{'data'}{$id}{'Cartn_z'} = $cartn_z;
            $atom_site{'data'}{$id}{'occupancy'} = $occupancy;
            $atom_site{'data'}{$id}{'B_iso_or_equiv_esd'} = $B_iso_or_equiv_esd;
            $atom_site{'data'}{$id}{'pdbx_formal_charge'} = $pdbx_formal_charge;
            $atom_site{'data'}{$id}{'pdbx_PDB_model_num'} = $pdbx_PDB_model_num;
        }
    }

    if( %{ $atom_site{'data'} } ) {
        $atom_site{'attributes'} = [
            'group_PDB', 'id', 'type_symbol', 'label_atom_id', 'label_alt_id',
            'label_comp_id', 'label_asym_id', 'label_seq_id', 'Cartn_x',
            'Cartn_y', 'Cartn_z', 'occupancy', 'B_iso_or_equiv_esd',
            'pdbx_formal_charge', 'pdbx_PDB_model_num'
        ];
        $atom_site{'is_loop'} = 1;
    }

    return \%atom_site;
}

1;
