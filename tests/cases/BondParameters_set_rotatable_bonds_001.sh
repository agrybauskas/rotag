#!/bin/bash

perl <<"END"
    use strict;
    use warnings;

    use Data::Dumper;

    use ForceField::Parameters;

    $Data::Dumper::Sortkeys = 1;
    $Data::Dumper::Indent = 1;

    my $PARAMETERS = Parameters->new();
    my $bond_parameters = BondParameters->new( $PARAMETERS );

    print Dumper $PARAMETERS;
END
