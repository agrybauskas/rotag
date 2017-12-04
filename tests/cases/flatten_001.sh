#!/bin/bash
cd "$(dirname "$0")"

matrix_file=../inputs/matrices_001.dat

../programs/flatten ${matrix_file}
