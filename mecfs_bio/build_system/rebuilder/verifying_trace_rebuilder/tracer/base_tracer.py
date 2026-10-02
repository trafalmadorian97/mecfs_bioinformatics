from abc import ABC, abstractmethod

from mecfs_bio.build_system.asset.base_asset import Asset


class Tracer(ABC):
    @abstractmethod
    def __call__(self, a: Asset) -> str:
        pass
