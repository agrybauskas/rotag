#!/bin/bash
cd "$(dirname "$0")"

matrix=../input/matrix_001.dat

../programs/transpose < ${matrix}
