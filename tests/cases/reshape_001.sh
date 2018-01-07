#!/bin/bash

export PERL5LIB=$(dirname "$0")/../../lib

element_list_file=$(dirname "$0")/../inputs/element-list-001.dat
dimensions="2,3 2,2"

$(dirname "$0")/../scripts/reshape "${dimensions}" ${element_list_file}
