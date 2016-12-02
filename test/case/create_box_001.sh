#!/bin/bash
cd "$(dirname "$0")"
../programs/create_box "label_atom_id CA,C,CB,OG" "Cartn_x,Cartn_y,Cartn_z" < ../input/serine_001.cif
