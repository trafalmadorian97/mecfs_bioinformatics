"""
Optional overrides for the gctb genome-wide fine-mapping MCMC chain (gctb --gwfm RC).

These map to gctb's basic MCMC command-line options (documented at
https://gctbhub.cloud.edu.au/software/gctb/). Every field is optional: a value of
None leaves gctb's built-in default in place, so the assembled command only carries
the flags that were explicitly overridden. A main use is a fast dummy run before a
full multi-hour job.
"""

from attrs import frozen


@frozen
class GctbMcmcOptions:
    """MCMC chain overrides for the gctb --gwfm RC step.

    chain_length is the total number of MCMC iterations; burn_in the number of
    leading iterations discarded; thin the sampling interval for recorded draws;
    out_freq how often gctb reports intermediate progress; seed the RNG seed for
    reproducible runs. Each is None to keep gctb's default. The invariants below
    make an unusable combination impossible to construct.
    """

    chain_length: int | None = None
    burn_in: int | None = None
    thin: int | None = None
    out_freq: int | None = None
    seed: int | None = None

    def __attrs_post_init__(self) -> None:
        for name, value in (
            ("chain_length", self.chain_length),
            ("burn_in", self.burn_in),
            ("thin", self.thin),
            ("out_freq", self.out_freq),
        ):
            assert value is None or value > 0, (
                f"{name} must be a positive number of iterations, got {value}"
            )
        if self.chain_length is not None and self.burn_in is not None:
            assert self.burn_in < self.chain_length, (
                f"burn_in ({self.burn_in}) must be less than "
                f"chain_length ({self.chain_length})"
            )


def render_mcmc_option_flags(options: GctbMcmcOptions) -> list[str]:
    """Render the set overrides as a flat list of gctb command-line tokens.
    """
    flags: list[str] = []
    if options.chain_length is not None:
        flags += ["--chain-length", str(options.chain_length)]
    if options.burn_in is not None:
        flags += ["--burn-in", str(options.burn_in)]
    if options.thin is not None:
        flags += ["--thin", str(options.thin)]
    if options.out_freq is not None:
        flags += ["--out-freq", str(options.out_freq)]
    if options.seed is not None:
        flags += ["--seed", str(options.seed)]
    return flags
