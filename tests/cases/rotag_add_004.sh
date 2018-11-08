#!/bin/bash

pdbx_file=$(dirname "$0")/../inputs/surrounded/serine-001.cif

rotag_add -H ${pdbx_file}
