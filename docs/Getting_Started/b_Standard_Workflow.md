# Standard Workflow

This page describes the standard workflow for contributing to the mecfs bioinformatics repository.



## Create a new branch

Use [checkout](https://git-scm.com/docs/git-checkout) to create a new branch to store your changes.  

```bash
git checkout -b <branchname>
```

## Make changes

Implement your changes.  Typically, this will involve:

- Defining a new `Task` object in the [assets][mecfs_bio.assets] directory.  This can often be accomplished using an existing `Task` class, but if no `Task` class meets your needs, you can implement a new one in the [task][mecfs_bio.build_system.task] directory.
- Adding an [analysis script][mecfs_bio.analysis] to materialize the new asset you defined.

If you added a new task class, you will also want to add a unit test to the `test_mecfs_bio directory` to verify that it works as intended on some dummy data.  To read about the principles of unit testing see Winters et al.[@winters2020software] (The chapter on testing can be found [here](https://abseil.io/resources/swe-book/html/ch12.html)).  

## Update documentation

If any documentation changes are needed, edit the corresponding files in the `docs` directory, then view the results by running the following command and accessing [http://localhost:8000](http://localhost:8000) in your browser.

```
pixi r mkdocs serve
```

To check for issues, add the `strict` flag, which will make the build abort if there are any warnings.

```
pixi r mkdocs serve --strict
```

## Run linters, formatters, and tests
The `invoke green` command runs linters, formatters, and tests. Run it with 

```bash
pixi r invoke green
```

It is advisable to run `invoke green` and fix any detected errors prior to submitting a PR.

##  Commit changes, create a PR

[Commit](https://git-scm.com/docs/git-commit) your changes and [push](https://git-scm.com/docs/git-push) them
```
git add mecfs_bio test_mecfs_bio; git commit -m '<your message>';git push --set-upstream origin <branchname>
```
Finally, use the github interface to create a PR from your changes, and ask a reviewer for a review.

