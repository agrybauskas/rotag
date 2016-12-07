#!/bin/bash
cd "$(dirname "$0")"

cif_file=../input/5svd_002.cif

../programs/obtain_atom_site < ${cif_file}
