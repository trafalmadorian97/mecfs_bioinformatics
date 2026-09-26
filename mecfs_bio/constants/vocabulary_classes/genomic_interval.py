from attrs import frozen


@frozen(slots=True)
class GenomicInterval:
    chrom: int
    start: int
    end: int
