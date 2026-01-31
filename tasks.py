from pathlib import Path

from invoke import task

NEW_UNIT_TEST_PATH = Path("test_mecfs_bio/unit")
SRC_PATH = Path("mecfs_bio")
DOCS_PATH = Path("docs")


# dev tasks
@task
def test(c):
    print("Running unit and integration tests with pytest...")
    cmd = f"pixi r python   -m pytest --typeguard-packages={SRC_PATH}  {NEW_UNIT_TEST_PATH}"
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
    c.run("pixi r  mypy .", pty=True)


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


@task(
    pre=[
        lintfix,
        format,
        spellcheck_docs,
        spellcheck_src,
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


@task(pre=[install_r_packages, green])
def init(c):
    """
    Initial repo setup
    """
    pass
