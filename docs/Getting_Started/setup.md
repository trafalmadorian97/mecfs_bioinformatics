# Initial Setup

We use Pixi are our package manager.  Setup [Pixi](https://pixi.sh/dev/) as described [here](https://pixi.sh/dev/installation/),
or using the bash command

```bash
curl -fsSL https://pixi.sh/install.sh | sh
```
source your `bashrc` to make sure `pixi` is on your path:
```bash
source ~/.bashrc
```


To verify everything is setup correctly, you can start by running the linters and unit tests:
```
pixi r invoke green
```