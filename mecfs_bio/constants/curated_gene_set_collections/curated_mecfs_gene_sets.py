"""
This file contains a selection of gene sets chosen because they relate to theories of the pathogensis
of ME/CFS that have been put forward.  The idea is to use them for a targeted gene set analysis.

"""

from attrs import frozen


@frozen
class MSigDBGeneSetSpec:
    exact_source: str | None
    standard_name: str
    systematic_name: str | None = None
    comment: str | None = None


CURATED_POTENTIAL_MECFS_GENE_SETS_AUTONOMIC = [
    # autonomic nervous system
    MSigDBGeneSetSpec(
        exact_source="GO:0048483",
        standard_name="GOBP_AUTONOMIC_NERVOUS_SYSTEM_DEVELOPMENT",
        comment="Justification: ME/CFS is hypothesized by some authors to be related to the autonomic nervous system.",
    ),
    MSigDBGeneSetSpec(
        exact_source="GO:0048484",  #  url: http://amigo.geneontology.org/amigo/term/GO:0048484
        standard_name="GOBP_ENTERIC_NERVOUS_SYSTEM_DEVELOPMENT",
        comment="Justification: ME/CFS is hypothesized by some authors to be related to the autonomic nervous system. "
        " The enteric nervous system is a branch of the autonomic nervous system",
    ),  #
    MSigDBGeneSetSpec(
        exact_source="GO:0048485",
        standard_name="GOBP_SYMPATHETIC_NERVOUS_SYSTEM_DEVELOPMENT",  #
        comment="Justification: ME/CFS is hypothesized by some authors to be related to the autonomic nervous system. "
        " The sympathetic nervous system is a branch of the autonomic nervous system",
    ),
    MSigDBGeneSetSpec(
        exact_source="GO:0048486",
        standard_name="GOBP_PARASYMPATHETIC_NERVOUS_SYSTEM_DEVELOPMENT",
        comment="Justification: ME/CFS is hypothesized by some authors to be related to the autonomic nervous system. "
        " The parasympathetic nervous system is a branch of the autonomic nervous system",
    ),
]

CURATED_POTENTIAL_MECFS_GENE_SETS_AUTONOMIC_SYNAPSE = [
    MSigDBGeneSetSpec(
        exact_source="GO:0071875",
        standard_name="GOBP_ADRENERGIC_RECEPTOR_SIGNALING_PATHWAY",
        comment="Justification: Adrenergic synapses are found in the peripheral sympathetic nervous system."
        "Some theorize ME/CFS is related to autonomic nervous system dysfunction",
    ),
    MSigDBGeneSetSpec(
        exact_source="GO:0007271",
        standard_name="GOBP_SYNAPTIC_TRANSMISSION_CHOLINERGIC",
        comment="Justification: Ganglionic synapses in the autonomic nervous system are largely use acetylcholine as "
        "their neurotransmitter. ",
    ),
]

CURATED_POTENTIAL_MECFS_GENE_SETS_BROAD_NEUROTRANSMISSION = [
    MSigDBGeneSetSpec(
        exact_source="hsa04080",
        systematic_name="M13380",
        standard_name="KEGG_NEUROACTIVE_LIGAND_RECEPTOR_INTERACTION",
        comment="This is a broad gene set capturing many aspects of neurotransmission. It thus serves to complement the more specific gene sets above",
    )
]


CURATED_POTENTIAL_MECFS_GENE_SETS_DRUG_RESPONSE = [
    MSigDBGeneSetSpec(
        exact_source="GO:0001975",
        standard_name="GOBP_RESPONSE_TO_AMPHETAMINE",
        comment="Justification: Anecdotally, some ME/CFS patients report some temporary benefits from stimulants.  "
        "Additionally, amphetamines are sometimes used as treatment for narcolepsy.  Some authors have theorized"
        "a connection between narcolepsy and ME/CFS",
    ),
    MSigDBGeneSetSpec(
        exact_source="GO:0031000",
        standard_name="GOBP_RESPONSE_TO_CAFFEINE",
        comment="Justification: Anecdotally, some ME/CFS patients report some temporary benefits from caffeine,"
        "[though long term, some report that caffeine actually makes them worse]",
    ),
]

CURATED_POTENTIAL_MECFS_GENE_SETS_VASCULAR_SMOOTH_MUSCLE = [
    MSigDBGeneSetSpec(
        exact_source="hsa04270",
        standard_name="KEGG_VASCULAR_SMOOTH_MUSCLE_CONTRACTION",
        systematic_name="M9387",
        comment="David Systrom and colleagues have observed low peak artereo-venous O2 difference in ME/CFS patients"
        "One explanation for this would be a peripheral shunting effect where a failure of peripheral microvasculature"
        "control prevents the appropriate direction of oxygenated blood to exercising muscle",
    ),
    MSigDBGeneSetSpec(
        exact_source="GO:0045906",
        systematic_name="M23756",
        standard_name="GOBP_NEGATIVE_REGULATION_OF_VASOCONSTRICTION",
        comment="Conceivable, peripheral shunting could be caused by a disruption of the vasoconstriction system",
    ),
]


CURATED_POTENTIAL_MECFS_GENE_SETS_NITRIC_OXIDE = [
    MSigDBGeneSetSpec(
        exact_source="GO:0007263",
        systematic_name="M15820",
        standard_name="GOBP_NITRIC_OXIDE_MEDIATED_SIGNAL_TRANSDUCTION",
        comment="This pathway captures the downstream effects of nitric oxide signalling.  Since nitric oxide plays a key role"
        "in the control of blood flow and vasodilation, inclusion of this pathway is justified by Systrom et al.'s findings",
    ),
    MSigDBGeneSetSpec(
        exact_source=None,
        systematic_name="M11650",
        standard_name="BIOCARTA_NOS1_PATHWAY",
        comment="Another nitric oxide pathway",
    ),
    MSigDBGeneSetSpec(
        exact_source="WP1995",
        systematic_name="M39677",
        standard_name="WP_EFFECTS_OF_NITRIC_OXIDE",
        comment="Another gene set for downstream effects of nitric oxide",
    ),
]

CURATED_POTENTIAL_MECFS_GENE_SETS_NEUROPATHY = [
    MSigDBGeneSetSpec(
        standard_name="GOCC_VOLTAGE_GATED_SODIUM_CHANNEL_COMPLEX",
        systematic_name="M17379",
        exact_source="GO:0001518",
        comment="Sodium channels are key to synaptic signaling.  If we think ME/CFS is a disease of synaptic connection,"
        "we might expect to find signal here",
    ),
    MSigDBGeneSetSpec(
        standard_name="GOBP_SENSORY_PERCEPTION_OF_PAIN",
        systematic_name="M12273",
        exact_source="GO:0019233",
        comment="Many ME/CFS patients have chronic pain.  If we believe ME/CFS is neurological, then perhaps the pain ME/CFS"
        "patients experience is neurological in origin",
    ),
    MSigDBGeneSetSpec(
        standard_name="GOBP_AXON_EXTENSION",
        systematic_name="M16466",
        exact_source="GO:0048675",
        comment="Failure of axon extension can cause neuropathy.  If we think ME/CFS is neurological, this is a potential cause",
    ),
    MSigDBGeneSetSpec(
        standard_name="GOBP_ANTEROGRADE_AXONAL_TRANSPORT",
        systematic_name="M29064",
        exact_source="GO:0008089",
        comment="Failure of axonal transport is another potential cause of neuropathy",
    ),
    MSigDBGeneSetSpec(
        standard_name="KEGG_NEUROTROPHIN_SIGNALING_PATHWAY",
        systematic_name="M16763",
        exact_source="hsa04722",
        comment="Failure of nerve growth signalling can cause neuropathy",
    ),
    MSigDBGeneSetSpec(
        standard_name="BIOCARTA_NGF_PATHWAY",
        systematic_name="M7860",
        exact_source=None,
        comment="Same motivation as above, but focusing on a narrower subset of the nerve growth signalling pathways.",
    ),
    MSigDBGeneSetSpec(
        standard_name="HP_DISTAL_PERIPHERAL_SENSORY_NEUROPATHY",
        systematic_name="M36986",
        exact_source="HP:0007067",
        comment="Small nerve fiber neuropathy has been associated with ME/CFS in some studies",
    ),
    MSigDBGeneSetSpec(
        standard_name="HP_SENSORY_AXONAL_NEUROPATHY",
        systematic_name="M36282",
        exact_source="HP:0003390",
        comment="Small fiber neuropathy has been associated with ME/CFS in some studies.  This focuses specifically on axonal disease",
    ),
]

CURATED_POTENTIAL_MECFS_GENE_SETS_MITOCHONDRIA = [
    MSigDBGeneSetSpec(
        standard_name="HALLMARK_OXIDATIVE_PHOSPHORYLATION",
        systematic_name="M5936",
        exact_source=None,
    ),
    MSigDBGeneSetSpec(
        standard_name="HP_MITOCHONDRIAL_RESPIRATORY_CHAIN_DEFECTS",
        exact_source="HP:0200125",
        systematic_name="M38933",
    ),
    MSigDBGeneSetSpec(
        standard_name="WP_ELECTRON_TRANSPORT_CHAIN_OXPHOS_SYSTEM_IN_MITOCHONDRIA",
        exact_source="WP111",
        systematic_name="M39417",
    ),
    MSigDBGeneSetSpec(
        standard_name="GOBP_MITOCHONDRIAL_ELECTRON_TRANSPORT_UBIQUINOL_TO_CYTOCHROME_C",
        exact_source="GO:0006122",
        systematic_name="M12278",
    ),
    MSigDBGeneSetSpec(
        standard_name="GOBP_CELLULAR_RESPONSE_TO_OXYGEN_CONTAINING_COMPOUND",
        exact_source="GO:1901701",
        systematic_name="M11609",
    ),
    MSigDBGeneSetSpec(
        standard_name="KEGG_ARACHIDONIC_ACID_METABOLISM",
        systematic_name="M5410",
        exact_source="hsa00590",
    ),
]

CURATED_POTENTIAL_MECFS_GENE_SETS_INFLAMMATION = [
    MSigDBGeneSetSpec(
        standard_name="REACTOME_THE_NLRP3_INFLAMMASOME",
        exact_source="R-HSA-844456",
        systematic_name="M1063",
    ),
    MSigDBGeneSetSpec(
        standard_name="GOBP_INFLAMMASOME_MEDIATED_SIGNALING_PATHWAY",
        exact_source="GO:0141084",
        systematic_name="M46883",
    ),
    MSigDBGeneSetSpec(
        standard_name="PID_TRAIL_PATHWAY",
        systematic_name="M79",
        exact_source=None,
        comment="David Systrom reported a trend to higher levels of TRAIL during exercise in ME/CFS patients",
    ),
]

CURATED_POTENTIAL_MECFS_GENE_SETS_RESPONSE_TO_VIRUS = []
