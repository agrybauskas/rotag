package AtomSite;

use strict;
use warnings;

# ------------------------- Constructor and destructor ------------------------ #

sub new
{
    my ( $class, $options ) = @_;
    my $self = {};

    return bless $self, $class;
}

sub destroy
{
    my ( $self ) = @_;
}

# --------------------------------- Methods ----------------------------------- #

1;
