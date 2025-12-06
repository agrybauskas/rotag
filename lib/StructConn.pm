package StructConn;

use strict;
use warnings;

use Exporter qw( import );
our @EXPORT_OK = qw( create_struct_conn );

use Carp;

use Version qw( $VERSION );

our $VERSION = $VERSION;

# ----------------------------------------------------------------------------- #

sub create_struct_conn
{
    my ( $atom1, $atom2, $type, $id ) = @_;
    my %struct_conn = (
        "$type$id" => {
            'id' => "$type$id",
            'conn_type_id' => "$type",
            'ptnr1_label_seq_id' => $atom1->{'label_seq_id'},
            'ptnr1_label_asym_id' => $atom1->{'label_asym_id'},
            'pdbx_ptnr1_label_alt_id' => $atom1->{'label_alt_id'},
            'ptnr1_auth_seq_id' => $atom1->{'auth_seq_id'},
            'ptnr1_auth_asym_id' => $atom1->{'auth_asym_id'},
            '[local]_ptnr1_pdbx_PDB_model_num' => $atom1->{'pdbx_PDB_model_num'},
            'ptnr1_label_atom_id' => $atom1->{'label_atom_id'},
            'ptnr2_label_seq_id' => $atom2->{'label_seq_id'},
            'ptnr2_label_asym_id' => $atom2->{'label_asym_id'},
            'pdbx_ptnr2_label_alt_id' => $atom2->{'label_alt_id'},
            'ptnr2_auth_seq_id' => $atom2->{'auth_seq_id'},
            'ptnr2_auth_asym_id' => $atom2->{'auth_asym_id'},
            '[local]_ptnr2_pdbx_PDB_model_num' => $atom2->{'pdbx_PDB_model_num'},
            'ptnr2_label_atom_id' => $atom2->{'label_atom_id'},
        }
    );
    return \%struct_conn;
}

1;
