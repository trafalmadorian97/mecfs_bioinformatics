#  Adding Packages

You may find that the analysis you would like to perform requires access to a software package that is not currently part of the `mecfs_bioinformatics` environment.

## Adding PyPi Packages

If your software package is available on PyPi, you can add it with 

```bash
pixi add --pypi <your_package_name>
```

This will update `pyproject.toml` and `pixi.lock`.  You can use `git add pyproject.toml pixi.lock; git commit -m <your_message>` to add the environment update to your branch. 


## Adding bioconda/ conda-forge packages

If your package is instead present on [bioconda](https://bioconda.github.io/) or [conda-forge](https://conda-forge.org/), you can add it with

```bash
pixi add <your_package_name>
```

Again, commit the resulting changes to `pyproject.toml` and `pixi.lock`.

Note that numerous Python packages are available on both PyPi and bioconda/ conda-forge.  In such cases, it is preferable to install them from PyPi.  Having Python packages installed from both sources can cause problems for Pixi's dependency resolver.