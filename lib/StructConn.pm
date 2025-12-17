package StructConn;

use strict;
use warnings;

use Exporter qw( import );
our @EXPORT_OK = qw( counter_struct_conn
                     create_struct_conn );

use Carp;

use Version qw( $VERSION );

our $VERSION = $VERSION;

# ----------------------------------------------------------------------------- #

sub create_struct_conn
{
    my ( $atom1, $atom2, $type, $id ) = @_;
    # TODO: need to find more robust solution for assigning certain unknown
    # values.
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

sub counter_struct_conn
{
    my ( $struct_conn ) = @_;
    my %counter = (
        'covale' => 1,
        'disulf' => 1,
        'hydrog' => 1,
        'metalc' => 1,
    );
    my ( $covale_last_id ) =
        sort { $b <=> $a } grep { s/^covale//g } keys %{ $struct_conn };
    my ( $disulf_last_id ) =
        sort { $b <=> $a } grep { s/^disulf//g } keys %{ $struct_conn };
    my ( $hydrog_last_id ) =
        sort { $b <=> $a } grep { s/^hydrog//g } keys %{ $struct_conn };
    my ( $metalc_last_id ) =
        sort { $b <=> $a } grep { s/^metalc//g } keys %{ $struct_conn };
    if( defined $covale_last_id ) {
        $counter{'covale'} = $covale_last_id + 1;
    }
    if( defined $disulf_last_id ) {
        $counter{'disulf'} = $disulf_last_id + 1;
    }
    if( defined $hydrog_last_id ) {
        $counter{'hydrog'} = $hydrog_last_id + 1;
    }
    if( defined $metalc_last_id ) {
        $counter{'metalc'} = $metalc_last_id + 1;
    }
    return \%counter;
}

1;
