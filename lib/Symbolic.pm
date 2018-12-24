package Symbolic;

sub new
{
    my ( $class, $args ) = @_;
    my $self = { 'is_evaluated' => $args->{'is_evaluated'},
                 'symbols' => $args->{'symbols'},
                 'matrix' => $args->{'matrix'} };
    return bless $self, $class;
}

sub evaluate
{
    my ( $self, $values ) = @_;

    for my $symbol ( @{ $self->{'symbols'} } ) {
        if( ! exists $values->{$symbol} ) {
            die "$symbol value is not passed";
        }
    }

    eval{
        $self->{'matrix'} =
            $self->{'matrix'}->( map { $values->{$_} } @{ $self->{'symbols'} } );
        $self->{'is_evaluated'} = 1;
    };
    if( @$ ) {
        die "$_";
    }

    return;
}

sub set_symbols
{
    my ( $self, $symbols ) = @_;
    $self->{'symbols'} = $symbols;
    return;
}

1;
