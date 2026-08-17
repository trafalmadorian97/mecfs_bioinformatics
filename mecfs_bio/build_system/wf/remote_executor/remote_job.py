from collections.abc import Mapping, Sequence
from pathlib import Path, PurePath

from attrs import frozen

from mecfs_bio.util.format_verify.docker_image import is_valid_docker_image_reference

# No `use_spot` / `instance_type` (see design spec). If a named-instance override is
# added later, refactor RemoteJob.resources to `RemoteResources | ExplicitInstance`
# rather than adding optional fields, keeping invalid states unrepresentable.


@frozen
class RemoteResources:
    memory_gb: int
    vcpus: int
    disk_gb: int
    region: str | None = None

    def __attrs_post_init__(self) -> None:
        assert self.memory_gb > 0 and self.vcpus > 0 and self.disk_gb > 0, (
            "RemoteResources fields must be positive"
        )


@frozen
class RemoteJob:
    image: str
    commands: Sequence[str]
    input_files: Mapping[Path, PurePath]
    s3_inputs: Mapping[str, PurePath]
    output_files: Sequence[PurePath]
    resources: RemoteResources

    def __attrs_post_init__(self) -> None:
        assert self.commands, "RemoteJob.commands must be non-empty"
        assert self.output_files, "RemoteJob.output_files must be non-empty"
        # s3_inputs keys are intentionally NOT validated as s3:// URIs here:
        # LocalDockerRemoteExecutor accepts a local directory path as an
        # s3_inputs key (its test seam), so the s3-only requirement is enforced
        # in the SkyPilot executor's build_sky_task instead.
        assert is_valid_docker_image_reference(self.image), (
            f"RemoteJob.image {self.image!r} is not a valid docker image reference"
        )
