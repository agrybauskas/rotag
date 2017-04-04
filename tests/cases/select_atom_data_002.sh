#!/bin/bash
cd "$(dirname "$0")"

cif_file=../inputs/5svd_002.cif
atom_specifier="id 1,2,3,4,5,6,7,8,9"
data_specifier="label_atom_id"

../programs/select_atom_data "${atom_specifier}" \
			     "${data_specifier}" \
			     < ${cif_file} \
| sort
