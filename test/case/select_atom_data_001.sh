#!/bin/bash
cd "$(dirname "$0")"
../programs/select_atom_data "label_atom_id CA,C,CB,OG" "Cartn_x,Cartn_y,Cartn_z" < ../input/serine_001.cif
