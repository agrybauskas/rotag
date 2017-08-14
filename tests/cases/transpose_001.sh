#!/bin/bash
cd "$(dirname "$0")"

matrix_file=../inputs/matrix_001.dat

../programs/transpose ${matrix_file}
