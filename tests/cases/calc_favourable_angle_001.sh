#!/bin/bash

export PERL5LIB=$(dirname "$0")/../../lib

pdbx_dump_file=$(dirname "$0")/../inputs/amino-acids/serine-rotation-only-001.dump
atom_id="152"
interaction_specifier="label_atom_id N,C,CA,O"
energy_cutoff_atom=0
checkable_angles="1.0*pi;1.1*pi;1.2*pi;1.3*pi;1.4*pi;1.5*pi;1.6*pi;1.7*pi;1.8*pi;1.9*pi;2.0*pi"

$(dirname "$0")/../scripts/calc_favourable_angle \
               ${atom_id} \
	       ${checkable_angles} \
               "${interaction_specifier}" \
	       ${energy_cutoff_atom} \
	       ${pdbx_dump_file}
