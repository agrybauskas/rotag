#!/bin/bash
cd "$(dirname "$0")"

symbols="chi"
matrices=../input/matrices_001.dat

../programs/dot_product ${symbols} < ${matrices}
