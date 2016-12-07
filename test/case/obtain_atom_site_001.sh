#!/bin/bash
cd "$(dirname "$0")"

cif_file=../input/serine_001.cif

../programs/obtain_atom_site < ${cif_file}
