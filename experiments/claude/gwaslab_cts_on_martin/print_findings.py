import pandas as pd

from experiments.claude.gwaslab_cts_on_martin.run import DATAFRAMES


def go():
    df = pd.read_csv("experiments/claude/gwaslab_cts_on_martin/results/SUMMARY_gtex_brain.csv", index_col="dataframe")
    df=df.loc[list(DATAFRAMES.keys()),:]
    result = df.to_markdown(floatfmt=".3e")
    with open("experiments/claude/gwaslab_cts_on_martin/results/summary_gtex_brain.md", "w") as f:
        f.write(result)


if __name__ == '__main__':
    go()