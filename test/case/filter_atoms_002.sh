#!/bin/bash
cd "$(dirname "$0")"
../programs/filter_atoms "label_atom_id CA,CB" < ../input/5svd_002.cif
