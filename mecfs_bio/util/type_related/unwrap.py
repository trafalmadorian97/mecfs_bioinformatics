def unwrap[T](x: T | None) -> T:
    if x is None:
        raise ValueError("Attempted to unwrap a None")
    return x
