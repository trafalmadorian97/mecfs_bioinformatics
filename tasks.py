from pathlib import Path

from invoke import task

NEW_UNIT_TEST_PATH = Path("test_mecfs_bio/unit")
SRC_PATH = Path("mecfs_bio")


# dev tasks
@task
def test(c):
    print("Running unit and integration tests with pytest")
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


@task(pre=[lintfix, format, typecheck, test])
def green(c):
    pass


## Library setup tasks


@task
def install_tabix_ubuntu(c):
    c.run("sudo apt-get install tabix")
