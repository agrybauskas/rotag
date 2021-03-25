package Interactions;

use strict;
use warnings;

use Digest::MD5;

# # Example.
# $ctx = Digest::MD5->new;
# $ctx->add($data);
# $ctx->addfile($file_handle);
# $digest = $ctx->digest;
# $digest = $ctx->hexdigest;
# $digest = $ctx->b64digest;

# ------------------------- Constructors/Destructors -------------------------- #

#
# Tracks overall interactions: their values, what is and not calculated, atom
# list and etc. The overall structure is in this form:
#
# { interactions => {
#     hex1 => {
#       atom_ids => [ atom_id1, atom_id2, ... ],
#       unique_residue_ids => [ unique_residue_key1, unique_residue_key2, ... ],
#       energy_func => \&energy_func,
#       energy_value => 1.0
#     },
#     hex2 => {
#       ...
#     },
#     ...
#   },
#   atom_ids => {
#     id1 => [hex1, hex2, ...],
#     id2 => [hex1],
#     ...
#   },
#   unique_reisdue_keys => {
#     unique_residue_key1 => [hex1, hex2, ...],
#     unique_residue_key2 => [hex1],
#     ...
#   },
# }
#
#
#
# Input:
#     ... - ...
# Output:
#     ... - ...
#

sub new
{
    my ( $class, $args ) = @_;
    my $self = {
        'interactions' => {},
        'atom_ids' => {},
        'unique_reisdue_keys' => {},
    };
    return bless $self, $class;
}

# ----------------------------- Setters/Getters ------------------------------- #

# --------------------------------- Methods ----------------------------------- #

sub add
{
    my ( $self, $atom_site, $options ) = @_;
    my $no_calculations = $options->{'no_calculations'};

    # MD5 hex creation for memoization purposes.
    my $md5hex = Digest::MD5->new;
    $md5hex->add( join ',', sort ( 1, 2, 3 ) );

    $self->{'interactions'}{$md5hex};
}

1;
