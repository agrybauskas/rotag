#!/bin/bash

pdbx_file=$(dirname "$0")/../inputs/empty.cif

rotag_mutate -m '1:SER,chi1=3.14' ${pdbx_file} 2>&1
