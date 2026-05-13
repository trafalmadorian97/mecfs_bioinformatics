from attrs import frozen


@frozen
class MSigDBGeneSetSpec:
    """
    A reference to an MSIGDB gene set
    exact_source, standard_name, systemic_name: column values for looking up a gene set.  Setting one of these to None means
    it will not be used to retrieve the gene set.

    Comment: documentation for human consumption.  Not used in actual look up
    """

    exact_source: str | None
    standard_name: str
    systematic_name: str | None = None
    comment: str | None = None
