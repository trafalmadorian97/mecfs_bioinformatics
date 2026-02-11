from pathlib import Path

from mecfs_bio.constants.gh_constants import GH_REPO_NAME
from mecfs_bio.util.github_commands.upload_download import release_from_dir, download_release_to_dir, does_release_exist


def go():
    release_from_dir(
        release_tag="dummy_release",
        dir_path=Path("dummy_dir"),
        title="dummy_release",
        draft=False,
        notes="testing",
        repo_name=GH_REPO_NAME
    )
    download_release_to_dir(
        release_tag="dummy_release",
        dir_path=Path("dummy_dir_12"),
        repo_name=GH_REPO_NAME,

    )


def go2():
    dummy_exists=does_release_exist(
        release_tag="dummy_release",
        repo_name=GH_REPO_NAME
    )

    dummy2_exists=does_release_exist(
        release_tag="dummy_release2",
        repo_name=GH_REPO_NAME
    )
    print(f"dummy exists:{dummy_exists}.  Dummy 2 exists: {dummy2_exists}")

if __name__ == "__main__":
    go()