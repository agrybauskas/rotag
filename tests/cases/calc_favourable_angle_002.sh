#!/bin/bash

export PERL5LIB=$(dirname "$0")/../../lib

pdbx_dump_file=$(dirname "$0")/../inputs/synthetic/xaa-rotation-only-002.dump
atom_id="3"
interaction_specifier="label_atom_id N,C,CA,O"
energy_cutoff_atom=0
checkable_angles="0.1*pi;0.2*pi;0.3*pi;0.4*pi;0.5*pi;0.6*pi;0.7*pi;0.8*pi;0.9*pi;1*pi;1.1*pi;1.2*pi;1.3*pi;1.4*pi;1.5*pi;1.6*pi;1.7*pi;1.8*pi;1.9*pi;2*pi"

$(dirname "$0")/../scripts/calc_favourable_angle \
               ${atom_id} \
	           ${checkable_angles} \
               "${interaction_specifier}" \
	           ${energy_cutoff_atom} \
	           ${pdbx_dump_file}
