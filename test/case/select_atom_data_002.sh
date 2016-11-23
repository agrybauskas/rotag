#!/bin/bash
cd "$(dirname "$0")"
../programs/select_atom_data "id 1,2,3,4,5,6,7,8,9" "label_atom_id" < ../input/5svd_002.cif
