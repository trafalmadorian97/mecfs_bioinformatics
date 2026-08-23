"""Tests for the order-agnostic allele key helpers."""

import narwhals as nw
import polars as pl

from mecfs_bio.build_system.task.ppp_database.allele_key import (
    unordered_allele_key,
    unordered_allele_key_narwhals,
)

_EA = "ea"
_NEA = "nea"
_KEY = "key"

# Both orientations of a variant, plus equal alleles and multi-character indels.
_ALLELE_PAIRS = [
    ("A", "T"),
    ("T", "A"),
    ("C", "G"),
    ("G", "C"),
    ("A", "A"),
    ("AC", "A"),
    ("A", "AC"),
]


def test_unordered_allele_key_is_orientation_agnostic():
    frame = pl.DataFrame({_EA: ["A", "T"], _NEA: ["T", "A"]})
    keys = frame.with_columns(unordered_allele_key(_EA, _NEA).alias(_KEY))[
        _KEY
    ].to_list()
    assert keys[0] == keys[1]


def test_narwhals_and_polars_allele_keys_agree():
    frame = pl.DataFrame(
        {_EA: [a for a, _ in _ALLELE_PAIRS], _NEA: [b for _, b in _ALLELE_PAIRS]}
    )
    polars_keys = frame.with_columns(unordered_allele_key(_EA, _NEA).alias(_KEY))[
        _KEY
    ].to_list()
    narwhals_keys = (
        nw.from_native(frame)
        .with_columns(unordered_allele_key_narwhals(_EA, _NEA).alias(_KEY))
        .to_native()[_KEY]
        .to_list()
    )
    assert polars_keys == narwhals_keys
