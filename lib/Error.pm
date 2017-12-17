package Error;

use strict;
use warnings;

# --------------------------------- Error Handling ---------------------------- #

#
# Object that can be created during the error and handeled properly.
#

sub new
{
    my $class = shift;
    my $self = { @_ };
    bless( $self, $class );
    return $self;
}

1;
