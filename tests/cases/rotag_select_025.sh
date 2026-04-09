#!/bin/bash

pdbx_file=$(dirname "$0")/../inputs/libraries/mg-with-sidechains-with-connections-library-001.cif

rotag_select -t 'atomname MG' --related-data --tags '_atom_site,_struct_conn,_[local]_rotamer_bond_parameter,_[local]_rotamer_energy' ${pdbx_file}
