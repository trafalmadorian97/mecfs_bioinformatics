from attrs import frozen


@frozen
class GenomicInterval:
    chrom: int
    start: int
    end: int
