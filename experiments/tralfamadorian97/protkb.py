from UniProtMapper import ProtKB
from UniProtMapper.uniprotkb_fields import reviewed, organism_name

protkb = ProtKB()

def go():
    # Find reviewed human proteins
    query = reviewed(True) & organism_name("human")
    fields= [
        "accession",
        "id",
        "gene_names",
        "protein_name",
        "organism_name",
        "organism_id",
        "go_id",
        "go_p",
        "go_c",
        "go_f",
        "cc_subcellular_location",
        "cc_function",
        "xref_ensembl",
        "xref_reactome"
    ]
    result = protkb.get(query, fields=fields)
    import pdb; pdb.set_trace()
    print("yo")

if __name__ == "__main__":
    go()