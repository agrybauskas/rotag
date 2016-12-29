#!/bin/bash
cd "$(dirname "$0")"

symbols="x,y,z"
matrices=../input/matrices_001.dat

../programs/dot_product ${symbols} < ${matrices}
