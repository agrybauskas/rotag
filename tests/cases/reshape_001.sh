#!/bin/bash
cd "$(dirname "$0")"

element_file=../inputs/element_list_001.dat
dimensions="2,3 2,2"

../programs/reshape "${dimensions}" ${element_file}
