"""
Abstract base class for Rebuilders.  See Andrey Mokhov, Neil Mitchell, and Simon Peyton Jones. Build systems à la carte.
"""

from abc import ABC, abstractmethod
from pathlib import Path

from mecfs_bio.build_system.asset.base_asset import Asset
from mecfs_bio.build_system.rebuilder.fetch.base_fetch import Fetch
from mecfs_bio.build_system.rebuilder.metadata_to_path.base_meta_to_path import (
    MetaToPath,
)
from mecfs_bio.build_system.task.base_task import Task
from mecfs_bio.build_system.wf.base_wf import WF


class Rebuilder[Info](ABC):
    """
    Key Operations:
    - Decide whether a given asset is up-to-date using information from Info.
    - If the asset is up-to-date, return it together with Info.
    - If the asset is not up-to-date, bring it up-to-date, update Info, and return the new values of both.
    """

    @abstractmethod
    def rebuild(
        self,
        task: Task,
        asset: Asset | None,
        fetch: Fetch,
        wf: WF,
        info: Info,
        meta_to_path: MetaToPath,
    ) -> tuple[Asset, Info]:
        pass

    @classmethod
    @abstractmethod
    def save_info(cls, info: Info, path: Path):
        pass
