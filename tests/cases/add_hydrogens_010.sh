#!/bin/bash
cd "$(dirname "$0")"

pdbx_file=../inputs/proline_010.cif

../programs/add_hydrogens ${pdbx_file}
