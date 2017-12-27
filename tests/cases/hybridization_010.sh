#!/bin/bash
cd "$(dirname "$0")"

pdbx_file=../inputs/proline_010.cif

../programs/hybridization ${pdbx_file}
