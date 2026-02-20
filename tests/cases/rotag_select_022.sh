#!/bin/bash

pdbx_file=$(dirname "$0")/../inputs/hetatoms/mg-with-sidechains-001.cif

rotag_select -t 'resname ASP' -s 'mainchain || atomname CB' -I 'resname MG' ${pdbx_file}
