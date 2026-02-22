from pathlib import Path

import yaml
from attrs import define, field

from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.meta.meta import Meta
from mecfs_bio.build_system.serialization.converter import CONVERTER_FOR_SERIALIZATION
from mecfs_bio.util.file_io.atomic_write import atomic_write_yaml

TraceRecord = tuple[str, list[tuple[AssetId, str]]]


@define
class VerifyingTraceInfo:
    """
    Stores persistent data relevant to the verifying trace rebuilder
    Note that this is mutable
    trace_store: A map associating asset ids with:
       - The traces of their corresponding assets
       -  A list of the asset ids of the dependencies of these assets, together with the traces of these dependencies
    trace_store is used to decide when an asset needs to be regenerated.  An asset is regenerated
    when its trace, or the traces of its dependencies have changed
    """

    trace_store: dict[AssetId, TraceRecord]
    must_rebuild: set[AssetId] = field(factory=set)

    @classmethod
    def empty(cls) -> "VerifyingTraceInfo":
        return cls(
            trace_store={},
            must_rebuild=set(),
        )

    def serialize(self, path: Path) -> None:
        conv = CONVERTER_FOR_SERIALIZATION
        path.parent.mkdir(parents=True, exist_ok=True)
        unstructured = conv.unstructure(self.trace_store)
        atomic_write_yaml(path=path, data=unstructured)

    @classmethod
    def deserialize(cls, path: Path) -> "VerifyingTraceInfo":
        conv = CONVERTER_FOR_SERIALIZATION
        with open(path) as infile:
            unstructured = yaml.load(infile, Loader=yaml.FullLoader)

        trace_store = conv.structure(unstructured, dict[AssetId, TraceRecord])

        return cls(trace_store=trace_store, must_rebuild=set())


def update_verifying_trace_info_in_place(
    verifying_trace_info: VerifyingTraceInfo,
    meta: Meta,
    new_value_trace: str,
    deps_traced: list[tuple[AssetId, str]],
) -> None:
    new_rebuild_set = set(verifying_trace_info.must_rebuild) - {meta.asset_id}
    new_trace_store = dict(verifying_trace_info.trace_store)
    new_trace_store[meta.asset_id] = new_value_trace, deps_traced
    verifying_trace_info.must_rebuild = new_rebuild_set
    verifying_trace_info.trace_store = new_trace_store
