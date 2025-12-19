# Standard Workflow

This page describes the standard workflow used to make a contribution to the mecft_bioinformatics repository.

## Clone repository

If you have not done previously, start by [cloning](https://git-scm.com/docs/git-clone) the repository.


```bash
git clone git@github.com:trafalmadorian97/mecfs_bioinformatics.git
```


## Create a new branch

The next step is to create a new branch to store your changes.  You can acomplish this with the [checkout](https://git-scm.com/docs/git-checkout) command

```bash
git checkout -b my_branch
```

## Make changes

The next step is to actually implement your changes.  Typically, this will involve.

- Implementing a new task object in the [assets][mecfs_bio.assets] directory.  This will often involve reusing existing Task classes, but no task class meets your needs, you can implement a new one in the [task][mecfs_bio.build_system.task] directory.
- Adding an [analysis script][mecfs_bio.analysis] to materialize the new asset yoy define.

If you added a new task class, you will also want to add a unit test to the test_mecfs_bio directory to verify that it works as intended on some dummy data.  For more on how to write unit tests See Winters et al.[@winters2020software] (The chapter on testing can be found [here](https://abseil.io/resources/swe-book/html/ch12.html)).  

## Run linters, formatters, and tests
The `green` command runs linters, formatters, and tests. This fixes the code match formatting conventions, and verifies that the changes have not broken any existing functionality as specified in the tests

```bash
pixi r invoke green
```


##  Commit changes, create a PR

Commit your changes and push them
```
git add mecfs_bio test_mecfs_bio; git commit -m '<your message>'; git push
```
Finally, use the github interface to create a PR from your changes, and ask for a review.

