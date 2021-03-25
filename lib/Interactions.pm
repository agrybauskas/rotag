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
#     md5hex1 => {
#       'atom_ids' => [ atom1_id, atom2_id, ... ],
#       'unique_residue_ids' => [ unique_residue1_id, unique_residue2_id, ... ],
#     },
#     md5hex2 => {
#       ...
#     },
#     ...
#   }
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
    my $self = { 'interactions' => {} };
    return bless $self, $class;
}

# ----------------------------- Setters/Getters ------------------------------- #

sub add
{

}

# --------------------------------- Methods ----------------------------------- #

1;
