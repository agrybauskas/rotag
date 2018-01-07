#!/bin/bash

export PERL5LIB=$(dirname "$0")/../../lib

x_distance="2.0" # In Angstrom.
y_distance="2.0" # In Angstrom.
z_distance="2.0" # In Angstrom.

$(dirname "$0")/../scripts/translation ${x_distance} ${y_distance} ${z_distance}
