#!/bin/bash

pdbx_file=$(dirname "$0")/../inputs/1a2p-shortened.cif
data_identifier='_atom_site'

$(dirname "$0")/../scripts/obtain_pdbx_data ${data_identifier} ${pdbx_file}
