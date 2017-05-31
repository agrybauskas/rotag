#!/bin/bash
cd "$(dirname "$0")"

matrix=../inputs/matrix_001.dat

../programs/transpose ${matrix}
