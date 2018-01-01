#!/bin/bash
cd "$(dirname "$0")"

atom_names_file=../inputs/atom_names_001.dat

../programs/sort_by_priority ${atom_names_file}
