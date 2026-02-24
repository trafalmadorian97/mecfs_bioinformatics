from pathlib import Path

from invoke import task

NEW_UNIT_TEST_PATH = Path("test_mecfs_bio/unit")
SRC_PATH = Path("mecfs_bio")
DOCS_PATH = Path("docs")

USER_AGENT = '"Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36"'

PULL_FIGURE_SCRIPT_PATH = Path("mecfs_bio/figures/key_scripts/pull_figures.py")

FIGS_PATH = Path("docs/_figs/")
FIGS_PATTERN = "_figs"


# dev tasks
@task
def test(c):
    print("Running unit and integration tests with pytest...")
    cmd = f"pixi r python   -m pytest --typeguard-packages={SRC_PATH}  {NEW_UNIT_TEST_PATH}"
    print(cmd)
    c.run(cmd, pty=True)


@task
def test_debug(c):
    """
    Run tests with extra information printed.  Useful for debugging.
    """
    print("Running unit and integration tests with pytest...")
    cmd = f"pixi r python   -m pytest -s --typeguard-debug-instrumentation --typeguard-packages={SRC_PATH}  {NEW_UNIT_TEST_PATH}"
    print(cmd)
    c.run(cmd, pty=True)


@task
def format(c):
    """
    Format code
    """
    print("Formatting with ruff...")
    c.run("pixi r  ruff format", pty=True)


@task
def formatcheck(c):
    """
    Check for format errors
    """
    print("Checking format with black...")
    c.run("pixi r ruff format --check .")


@task
def lintfix(c):
    print("linting and applying lint auto-fixes using ruff...")
    c.run(" pixi r  ruff check --fix --unsafe-fixes")


@task
def lintcheck(c):
    print("linting using ruff...")
    c.run("pixi r ruff check")


@task
def typecheck(c):
    """
    Check for type errors
    """
    print("Typechecking with mypy...")
    c.run(f"pixi r  mypy {SRC_PATH} {NEW_UNIT_TEST_PATH}", pty=True)


@task
def spellcheck_docs(c):
    """
    check docs for spelling errors
    Edit _typos.toml to add exceptions
    """
    print(
        f"Checking documentation spelling using typos... (Add exceptions to {DOCS_PATH}/_typos.toml)"
    )
    c.run(f"pixi r  typos {DOCS_PATH}", pty=True)


@task
def spellcheck_src(c):
    """
    check code for spelling errors
    Edit _typos.toml to add exceptions
    """
    print(
        f"Checking source spelling using typos... (Add exceptions to {SRC_PATH}/_typos.toml)"
    )
    c.run(f"pixi r  typos {SRC_PATH}", pty=True)


@task
def checkimports(c):
    """
    Use import linter to enforce architectural constraints
    """
    print("Checking architectural constraints using import-linter...")
    c.run("pixi r lint-imports", pty=True)


@task
def check_all_links(c):
    """
    Check all links with lychee

    I added "403" to the list of acceptable status codes to prevent this check from failing due to anti-bot systems.
    """
    print("Checking links with lychee...")
    c.run(
        f"pixi r lychee --insecure --cache --accept 100..=103,200..=299,403  --cache-exclude-status 400..=999 --exclude {FIGS_PATTERN} --user-agent {USER_AGENT}  {SRC_PATH} {DOCS_PATH}  "
    )


@task
def check_local_links(c):
    """
    Check local links with lychee
    """
    print("Checking offline links with lychee...")
    cmd = f"pixi r lychee --exclude {FIGS_PATTERN}  --offline {SRC_PATH} {DOCS_PATH}"
    print(f"running {cmd}")
    c.run(cmd)


@task(
    pre=[
        lintfix,
        format,
        spellcheck_docs,
        spellcheck_src,
        check_local_links,
        checkimports,
        typecheck,
        test,
    ]
)
def green(c):
    pass


@task
def install_r_packages(c):
    """
    Install R packages such as TwoSampleMR
    """
    c.run("pixi r install-mr", pty=True)


### Figures and Documentation


@task()
def pfig(c):
    """
    Pull figures from github to populate the _figs directory
    """
    c.run(f"pixi r python {PULL_FIGURE_SCRIPT_PATH}")


@task
def serve_docs(c, strict: bool = True):
    """
    Use mkdocs to serve documentation
    """
    cmd = "pixi r mkdocs serve"
    if strict:
        cmd += " --strict"
    print("Serving documentation...")
    print(f"running {cmd}")
    c.run(cmd, pty=True)


@task(pre=[pfig, serve_docs])
def sdocs(c):
    """
    Retrieve figures, then serve docs
    """


# initialization
@task(pre=[install_r_packages, pfig, green])
def init(c):
    """
    Initial repo setup
    """
    pass
