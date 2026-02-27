import pandas as pd


def make_bed_file():
    df = pd.read_csv("assets/base_asset_store/reference_data/magma_reference_data/ensembl_gene_locations/raw/ENSGv102.coding.genes.txt",sep="\s+",
                     header=None, names=["ensembl", "chrom", "start", "end","strand", "name"])
    df["intensity"]=100
    df.loc[:,["chrom","start","end","name", "intensity","strand"]].to_csv(
        "data/dummy_gene_track_file.txt",sep="\t",index=False,
        header=False
    )


if __name__ == "__main__":
    make_bed_file()
