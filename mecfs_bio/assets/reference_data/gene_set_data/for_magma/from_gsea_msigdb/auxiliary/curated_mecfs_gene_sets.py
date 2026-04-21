from attrs import frozen


@frozen
class MSigDBGeneSetSpec:
   exact_source:str|None
   standard_name: str
   systematic_name:str|None=None
   comment: str | None=None


CURATED_POTENTIAL_MECFS_GENE_SETS_AUTONOMIC = [

    # autonomic nervous system
    MSigDBGeneSetSpec(
        exact_source="GO:0048483",
        standard_name="GOBP_AUTONOMIC_NERVOUS_SYSTEM_DEVELOPMENT",
        comment="Justification: ME/CFS is hypothesized by some authors to be related to the autonomic nervous system."
    ),
    MSigDBGeneSetSpec(
    exact_source="GO:0048484",#  url: http://amigo.geneontology.org/amigo/term/GO:0048484
    standard_name="GOBP_ENTERIC_NERVOUS_SYSTEM_DEVELOPMENT",
        comment="Justification: ME/CFS is hypothesized by some authors to be related to the autonomic nervous system. "
                " The enteric nervous system is a branch of the autonomic nervous system"
    ), #
    MSigDBGeneSetSpec(
        exact_source="GO:0048485",
        standard_name="GOBP_SYMPATHETIC_NERVOUS_SYSTEM_DEVELOPMENT", #
        comment="Justification: ME/CFS is hypothesized by some authors to be related to the autonomic nervous system. "
                " The sympathetic nervous system is a branch of the autonomic nervous system"
    ),
    MSigDBGeneSetSpec(
            exact_source="GO:0048486",
        standard_name="GOBP_PARASYMPATHETIC_NERVOUS_SYSTEM_DEVELOPMENT",
        comment="Justification: ME/CFS is hypothesized by some authors to be related to the autonomic nervous system. "
                " The parasympathetic nervous system is a branch of the autonomic nervous system"
    ),
]

CURATED_POTENTIAL_MECFS_GENE_SETS_AUTONOMIC_SYNAPSE = [
    MSigDBGeneSetSpec(
        exact_source="GO:0071875",
        standard_name="GOBP_ADRENERGIC_RECEPTOR_SIGNALING_PATHWAY",
        comment="Justification: Adrenergic synapses are found in the peripheral sympathetic nervous system."
                "Some theorize ME/CFS is related to autonomic nervous system dysfunction"
    ),

    MSigDBGeneSetSpec(
        exact_source="GO:0007271",
        standard_name="GOBP_SYNAPTIC_TRANSMISSION_CHOLINERGIC",
        comment="Justification: Ganglionic synapses in the autonomic nervous system are largely use acetylcholine as "
                "their neurotransmitter. "
    ),

]


CURATED_POTENTIAL_MECFS_GENE_SETS_DRUG_RESPONSE=[

    MSigDBGeneSetSpec(
        exact_source="GO:0001975",
        standard_name="GOBP_RESPONSE_TO_AMPHETAMINE",
        comment="Justification: Anecdotally, some ME/CFS patients report some temporary benefits from stimulants.  "
                "Additionally, amphetamines are sometimes used as treatment for narcolepsy.  Some authors have theorized"
                "a connection between narcolepsy and ME/CFS"
    ),

    MSigDBGeneSetSpec(
        exact_source="GO:0031000",
        standard_name="GOBP_RESPONSE_TO_CAFFEINE",
        comment="Justification: Anecdotally, some ME/CFS patients report some temporary benefits from caffeine,"
                "[though long term, some report that caffeine actually makes them worse]"

    ),

]

CURATED_POTENTIAL_MECFS_GENE_SETS_VASCULAR_SMOOTH_MUSCLE=[]

CURRATED_POTENTIAL_MECFS_GENE_SETS_RESPONSE_TO_VIRUS=[]