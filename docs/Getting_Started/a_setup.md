# Initial Setup


If you have not done so previously, [clone](https://git-scm.com/docs/git-clone) the repository.


```bash
git clone git@github.com:trafalmadorian97/mecfs_bioinformatics.git
```

We use [Pixi](https://pixi.sh/dev/) as are our package manager.  Setup Pixi as described [here](https://pixi.sh/dev/installation/),
using the bash command

```bash
curl -fsSL https://pixi.sh/install.sh | sh
```

After setup, source your `bashrc` to make sure `pixi` is on your path:

```bash
source ~/.bashrc
```

To verify everything is working correctly, try: 

```bash
pixi r invoke init
```

This will install additional required packages, then run linters and tests.


A natural next step is to run the basic analysis of the Decode ME data using the script [here][mecfs_bio.analysis.decode_me_initial_analysis].

You can use pixi to run this or another script via the command:

```bash
pixi r python [PATH_TO_SCRIPT]
```