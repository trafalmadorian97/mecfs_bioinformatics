"""
Named outcomes of genome-reference harmonization: allele actions, palindrome decisions
and drop reasons.

Each Literal alias is what the type checker sees; the constants are what code uses, so
no outcome string is typed more than twice (once in the alias, once in its constant).
An import-time assertion keeps each alias and its constants in step. A StrEnum cannot
replace them, because polars converts Python Enum members into its own Enum dtype
inside expressions.
"""

from typing import Final, Literal, get_args

AlleleAction = Literal["keep", "swap", "complement", "complement_swap"]
ACTION_KEEP: Final[AlleleAction] = "keep"
ACTION_SWAP: Final[AlleleAction] = "swap"
ACTION_COMPLEMENT: Final[AlleleAction] = "complement"
ACTION_COMPLEMENT_SWAP: Final[AlleleAction] = "complement_swap"
assert set(get_args(AlleleAction)) == {
    ACTION_KEEP,
    ACTION_SWAP,
    ACTION_COMPLEMENT,
    ACTION_COMPLEMENT_SWAP,
}

PalindromeDecision = Literal["keep", "strand_flip", "unresolved"]
PALINDROME_KEEP: Final[PalindromeDecision] = "keep"
PALINDROME_STRAND_FLIP: Final[PalindromeDecision] = "strand_flip"
PALINDROME_UNRESOLVED: Final[PalindromeDecision] = "unresolved"
assert set(get_args(PalindromeDecision)) == {
    PALINDROME_KEEP,
    PALINDROME_STRAND_FLIP,
    PALINDROME_UNRESOLVED,
}

DropReason = Literal[
    "invalid_allele",
    "not_on_reference",
    "indel_not_on_reference",
    "palindromic_mnp_untrusted",
    "palindrome_unresolved",
    "ambiguous_indel_no_eaf",
    "ambiguous_indel_not_in_panel",
    "ambiguous_indel_af_mismatch",
    "ambiguous_indel_af_indecisive",
]
DROP_INVALID_ALLELE: Final[DropReason] = "invalid_allele"
DROP_NOT_ON_REFERENCE: Final[DropReason] = "not_on_reference"
DROP_INDEL_NOT_ON_REFERENCE: Final[DropReason] = "indel_not_on_reference"
DROP_PALINDROMIC_MNP_UNTRUSTED: Final[DropReason] = "palindromic_mnp_untrusted"
DROP_PALINDROME_UNRESOLVED: Final[DropReason] = "palindrome_unresolved"
DROP_AMBIGUOUS_INDEL_NO_EAF: Final[DropReason] = "ambiguous_indel_no_eaf"
DROP_AMBIGUOUS_INDEL_NOT_IN_PANEL: Final[DropReason] = "ambiguous_indel_not_in_panel"
DROP_AMBIGUOUS_INDEL_AF_MISMATCH: Final[DropReason] = "ambiguous_indel_af_mismatch"
DROP_AMBIGUOUS_INDEL_AF_INDECISIVE: Final[DropReason] = "ambiguous_indel_af_indecisive"
assert set(get_args(DropReason)) == {
    DROP_INVALID_ALLELE,
    DROP_NOT_ON_REFERENCE,
    DROP_INDEL_NOT_ON_REFERENCE,
    DROP_PALINDROMIC_MNP_UNTRUSTED,
    DROP_PALINDROME_UNRESOLVED,
    DROP_AMBIGUOUS_INDEL_NO_EAF,
    DROP_AMBIGUOUS_INDEL_NOT_IN_PANEL,
    DROP_AMBIGUOUS_INDEL_AF_MISMATCH,
    DROP_AMBIGUOUS_INDEL_AF_INDECISIVE,
}
