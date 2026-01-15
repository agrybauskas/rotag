#!/bin/bash

pdbx_file=$(dirname "$0")/../inputs/hetatoms/1eo3-active-site-with-hetatoms-selected-001.cif

rotag_connect --metalc '1,2' ${pdbx_file}
