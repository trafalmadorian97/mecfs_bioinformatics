import pytest

from mecfs_bio.build_system.task.sbayesrc.gctb_mcmc_options import (
    GctbMcmcOptions,
    render_mcmc_option_flags,
)


def test_all_default_options_render_to_no_flags() -> None:
    assert render_mcmc_option_flags(GctbMcmcOptions()) == []


def test_set_options_render_as_flag_value_token_pairs() -> None:
    flags = render_mcmc_option_flags(
        GctbMcmcOptions(chain_length=100, burn_in=20, thin=5, out_freq=10, seed=123)
    )
    # Rendered as a flat token list ready to splice into the command line.
    assert flags == [
        "--chain-length",
        "100",
        "--burn-in",
        "20",
        "--thin",
        "5",
        "--out-freq",
        "10",
        "--seed",
        "123",
    ]


def test_burn_in_must_be_less_than_chain_length() -> None:
    with pytest.raises(AssertionError):
        GctbMcmcOptions(chain_length=100, burn_in=100)


def test_iteration_counts_must_be_positive() -> None:
    with pytest.raises(AssertionError):
        GctbMcmcOptions(chain_length=0)
