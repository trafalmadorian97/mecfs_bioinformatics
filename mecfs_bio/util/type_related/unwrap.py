from typing import Sequence


def unwrap[T](x: T | None) -> T:
    if x is None:
        raise ValueError("Attempted to unwrap a None")
    return x


def mapped_unwrap[T](x: Sequence[T | None]) -> Sequence[T]:
    return [unwrap(y) for y in x]
