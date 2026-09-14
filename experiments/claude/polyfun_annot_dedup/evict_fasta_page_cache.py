"""Evict ~/.gwaslab/hg19.fa from the OS page cache (no root needed) so the next benchmark is cold."""
import os
from pathlib import Path

fd = os.open(Path.home() / ".gwaslab/hg19.fa", os.O_RDONLY)
os.posix_fadvise(fd, 0, 0, os.POSIX_FADV_DONTNEED)
os.close(fd)
print("evicted hg19.fa from page cache")
