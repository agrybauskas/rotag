#!/bin/bash
cd "$(dirname "$0")"

matrices=../input/matrices_001.dat

../programs/dot_product < ${matrices}
