#!/bin/bash

export PERL5LIB=$(dirname "$0")/../../lib

angle_range="0-2,4-5"
small_angle=0.2

$(dirname "$0")/../scripts/sample_angles "${angle_range}" "${small_angle}"
