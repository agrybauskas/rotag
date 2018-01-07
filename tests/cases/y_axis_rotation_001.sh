#!/bin/bash

export PERL5LIB=$(dirname "$0")/../../lib

angle="1.5708" # In radians.

$(dirname "$0")/../scripts/y_axis_rotation ${angle}
