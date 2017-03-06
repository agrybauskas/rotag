#!/bin/bash
cd "$(dirname "$0")"

cif_file=../inputs/serine_001.cif
atom_specifier="label_seq_id 18"

../programs/cf_serine "${atom_specifier}" < ${cif_file}
