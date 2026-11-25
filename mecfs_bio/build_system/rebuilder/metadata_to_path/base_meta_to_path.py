from abc import ABC, abstractmethod
from pathlib import Path

from mecfs_bio.build_system.meta.meta import Meta


class MetaToPath(ABC):
    @abstractmethod
    def __call__(self, m: Meta) -> Path:
        pass
