#!/bin/bash
cd "$(dirname "$0")"

pdbx_file=../inputs/aspartic_acid_006.cif
potential="leonard_jones"
cutoff=0.04
target_resi="219"
small_angle="0.4*pi"

../programs/draw_density "${potential}" "${cutoff}" "${target_resi}" \
			 "${small_angle}" ${pdbx_file}
