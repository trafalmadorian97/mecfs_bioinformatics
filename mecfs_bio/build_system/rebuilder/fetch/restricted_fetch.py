from attrs import frozen

from mecfs_bio.build_system.asset.base_asset import Asset
from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.rebuilder.fetch.base_fetch import Fetch
from mecfs_bio.build_system.task.base_task import Task


@frozen(slots=True)
class RestrictedFetch(Fetch):
    """
    A Fetch the prevents the retrieval of assets except those on a given list.
    Used to ensure the declared dependencies of a task are accurate.
    """

    inner: Fetch
    meta_deps: frozenset[AssetId]

    def __call__(self, asset_id: AssetId) -> Asset:
        if asset_id not in self.meta_deps:
            raise ValueError(
                f"Attempted to fetch asset {asset_id}, but only assets {self.meta_deps} are declared as dependencies."
            )
        return self.inner(asset_id)

    @classmethod
    def from_task(cls, fetch: Fetch, task: Task):
        return cls(fetch, frozenset({dep.meta.asset_id for dep in task.deps}))
