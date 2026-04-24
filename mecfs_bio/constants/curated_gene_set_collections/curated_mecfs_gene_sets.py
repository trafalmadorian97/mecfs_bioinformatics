"""
This file contains a selection of gene sets chosen because they relate to theories of the pathogensis
of ME/CFS that have been put forward.  The idea is to use them for a targeted gene set analysis.

"""

from mecfs_bio.constants.vocabulary_classes.gene_set import MSigDBGeneSetSpec

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
        "The enteric nervous system is a branch of the autonomic nervous system",
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
        exact_source="GO:0019229",
        systematic_name="M14342",
        standard_name="GOBP_REGULATION_OF_VASOCONSTRICTION",
        comment="Peripheral shunting could be caused by a disruption of the vasoconstriction system",
    ),
    MSigDBGeneSetSpec(
        exact_source="GO:0008217",
        systematic_name="M10469",
        standard_name="GOBP_REGULATION_OF_BLOOD_PRESSURE",
        comment="constriction and dilation of blood vessels is key to the control of blood pressure",
    ),
]

CURATED_POTENTIAL_MECFS_GENE_SETS_ENDOTHELIAL = [
    MSigDBGeneSetSpec(
        standard_name="GOBP_ENDOTHELIAL_CELL_ACTIVATION",
        exact_source="GO:0042118",
        systematic_name="M23506",
        comment="The Wirth-Scheibenbogen theory postulates that endothelial damage is a key driver",
    )
]

CURATED_POTENTIAL_MECFS_GENE_SETS_ANGIOGENSIS = [
    MSigDBGeneSetSpec(
        standard_name="HALLMARK_ANGIOGENESIS",
        systematic_name="M5944",
        exact_source=None,
        comment="Some long covid researchers report elevation of angiogensis factors in early stages of long covid",
    ),
    MSigDBGeneSetSpec(
        standard_name="HALLMARK_HYPOXIA",
        systematic_name="M5891",
        exact_source=None,
        comment="hypoxia is a primary driver of angiogensis",
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
        comment="Several authors have postulated a role for mitochondrial dysfunction in ME/CFS, though the"
        "evidence for this has been patchy.  Oxidative phosphorylation is the major energy producing action of mitochondria.",
    ),
    MSigDBGeneSetSpec(
        standard_name="HP_MITOCHONDRIAL_RESPIRATORY_CHAIN_DEFECTS",
        exact_source="HP:0200125",
        systematic_name="M38933",
        comment="see above",
    ),
    MSigDBGeneSetSpec(
        standard_name="GOBP_CELLULAR_RESPONSE_TO_OXYGEN_CONTAINING_COMPOUND",
        exact_source="GO:1901701",
        systematic_name="M11609",
        comment="see above",
    ),
    MSigDBGeneSetSpec(
        standard_name="KEGG_ARACHIDONIC_ACID_METABOLISM",
        systematic_name="M5410",
        exact_source="hsa00590",
        comment="see above",
    ),
    MSigDBGeneSetSpec(
        standard_name="HALLMARK_GLYCOLYSIS",
        systematic_name="M5937",
        exact_source=None,
        comment="Glycolysis is the less-efficient but anaerobic alternative to the glycolysis pathway",
    ),
    MSigDBGeneSetSpec(
        standard_name="KEGG_CITRATE_CYCLE_TCA_CYCLE",
        systematic_name="M3985",
        exact_source="hsa00020",
        comment="The citric acid cycle is an immediate precursor step to oxidative phosphorylation",
    ),
]

CURATED_POTENTIAL_MECFS_GENE_SETS_NEUROINFLAMMATION = [
    MSigDBGeneSetSpec(
        standard_name="GOBP_MICROGLIAL_CELL_ACTIVATION_INVOLVED_IN_IMMUNE_RESPONSE",
        systematic_name="M45050",
        exact_source="GO:0002282",
        comment="Researchers like Jarred Younger have hypothesis that some form of immune activation in the central nervous "
        "system could drive ME/CFS",
    ),
    MSigDBGeneSetSpec(
        standard_name="GOBP_ASTROCYTE_ACTIVATION",
        systematic_name="M23832",
        exact_source="GO:0048143",
        comment="Researchers like Jarred Younger have hypothesis that some form of immune activation in the central nervous "
        "system could drive ME/CFS",
    ),
]

CURATED_POTENTIAL_MECFS_GENE_SETS_INFLAMMATION = [
    MSigDBGeneSetSpec(
        standard_name="REACTOME_THE_NLRP3_INFLAMMASOME",
        exact_source="R-HSA-844456",
        systematic_name="M1063",
    ),
    MSigDBGeneSetSpec(
        standard_name="PID_TRAIL_PATHWAY",
        systematic_name="M79",
        exact_source=None,
        comment="David Systrom reported a trend to higher levels of TRAIL during exercise in ME/CFS patients",
    ),
    MSigDBGeneSetSpec(
        standard_name="HALLMARK_COMPLEMENT",
        systematic_name="M5921",
        exact_source=None,
        comment="Akiko Iwasaki argues for a role of the complement system in long covid.",
    ),
]

CURATED_POTENTIAL_MECFS_GENE_SETS_ADAPTIVE_IMMUNITY = [
    MSigDBGeneSetSpec(
        standard_name="GOBP_ANTIGEN_PROCESSING_AND_PRESENTATION",
        systematic_name="M12794",
        exact_source="GO:0019882",
        comment="Several theories classify ME/CFS as an unusual form of autoimmune disease."
        " This broad gene set captures captures the process wby which antigens are prepared are presented to T cells,"
        "which is altered in some autoimmune conditions.",
    ),
    MSigDBGeneSetSpec(
        standard_name="REACTOME_MHC_CLASS_II_ANTIGEN_PRESENTATION",
        systematic_name="M705",
        exact_source="R-HSA-2132295",
        comment="This is a more granular antigen presentation gene set focused on presentation by MHC 2, which is found"
        "primarily on antigen-presenting cells",
    ),
    MSigDBGeneSetSpec(
        standard_name="REACTOME_TCR_SIGNALING",
        systematic_name="M15381",
        exact_source="R-HSA-202403",
        comment="A finer-grained autoimmune gene set, focusing on the T cell axis.  Some authors have noted that "
        "post-infectious autoimmune conditions like reactive arthritis are often T-cell driven.  Thus, given that "
        "ME/CFS is seen as post-infectious, it may make sense to look for T cell signal",
    ),
    MSigDBGeneSetSpec(
        standard_name="REACTOME_SIGNALING_BY_THE_B_CELL_RECEPTOR_BCR",
        systematic_name="M608",
        exact_source="R-HSA-983705",
        comment="Alternatively, a group of researchers has theorized that autoantibodies to G-Protein-Coupled receptors"
        "may be central to ME/CFS pathogeneisis.  If this is correct, we might expect to see signals related to B-cell activation,"
        " as are found in this B-cell receptor gene set",
    ),
    MSigDBGeneSetSpec(
        standard_name="REACTOME_REGULATION_OF_T_CELL_ACTIVATION_BY_CD28_FAMILY",
        systematic_name="M17386",
        exact_source="R-HSA-388841",
        comment="T cells require co-stimulatory signals in addition to activation of their TCR. "
        "This requirement for a costimulatory check that maintains self tolerance.  Failure on this pathway "
        "can lead to autoimmunity",
    ),
    MSigDBGeneSetSpec(
        standard_name="WP_TH17_CELL_DIFFERENTIATION_PATHWAY",
        systematic_name="M45541",
        exact_source="WP5130",
        comment="Thss is a narrower pathway that is important so some autoimmune diseases (seronegative spondylo-arthritis) but not others"
        "(e.g. rheumatoid arthritis).  If do find autoimmune signal, it will be interested to see on which side of the divide it lies",
    ),
]

CURATED_POTENTIAL_MECFS_GENE_SETS_RESPONSE_TO_VIRUS = [
    MSigDBGeneSetSpec(
        standard_name="GOBP_RESPONSE_TO_VIRUS",
        systematic_name="M16779",
        exact_source="GO:0009615",
        comment="ME/CS is often described as a post-viral condition, and many patients report an onset after a severe infection.",
    ),
    MSigDBGeneSetSpec(
        standard_name="HALLMARK_INTERFERON_GAMMA_RESPONSE",
        systematic_name="M5913",
        exact_source=None,
        comment="Certain authors have observed that ME/CFS is often post-viral, that interferon gamma is part of the anti-viral response, "
        "And that interferon gamma induces severe fatigue",
    ),
    MSigDBGeneSetSpec(
        standard_name="HALLMARK_INTERFERON_ALPHA_RESPONSE",
        systematic_name="M5911",
        exact_source=None,
        comment="Interferon alpha is another branch of the anti-viral response",
    ),
]
CURATED_POTENTIAL_MECFS_GENE_SETS_HPA_AXIS = [
    MSigDBGeneSetSpec(
        standard_name="WP_GLUCOCORTICOID_RECEPTOR_PATHWAY",
        systematic_name="M40042",
        exact_source="WP2880",
        comment="Some author argue MECFS involves a blunted HPA axis",
    ),
    MSigDBGeneSetSpec(
        standard_name="WP_CORTICOTROPINRELEASING_HORMONE_SIGNALING",
        systematic_name="M39441",
        exact_source="WP2355",
        comment="A Dutch autopsy study of the brains of ME/CFS patients reports reported highly depleted "
        "CRH neurons",
    ),
    MSigDBGeneSetSpec(
        standard_name="GOBP_POSITIVE_REGULATION_OF_CORTICOSTEROID_HORMONE_SECRETION",
        systematic_name="M25536",
        exact_source="GO:2000848",
    ),
]

CURATED_POTENTIAL_MECFS_GENE_SETS_CLOTTING = [
    MSigDBGeneSetSpec(
        standard_name="HALLMARK_COAGULATION",
        systematic_name="M5946",
        exact_source=None,
        comment="Pretorius argues for an important role for microclots",
    )
]


CURATED_POTENTIAL_MECFS_GENE_SETS_METABOLIC_TRAP = [
    MSigDBGeneSetSpec(
        standard_name="KEGG_TRYPTOPHAN_METABOLISM",
        systematic_name="M980",
        exact_source="hsa00380",
        comment="I don't understand the details of the metabolic trap hypothesis, but according to a cursory search,"
        "it relates to tryptophan metabolism",
    ),
]

CURATED_POTENTIAL_MECFS_GENE_SETS_MAST_CELL_ACTIVATION = [
    MSigDBGeneSetSpec(
        standard_name="GOBP_MAST_CELL_ACTIVATION",
        systematic_name="M12306",
        exact_source="GO:0045576",
        comment="There is a long history of speculation that ME/CFS is linked to mast cell activation, though the evidence for this "
        "is generally low quality",
    )
]

CURATED_POTENTIAL_MECFS_GENE_SETS_SLEEP_DISTURBANCE = [
    MSigDBGeneSetSpec(
        standard_name="KEGG_CIRCADIAN_RHYTHM_MAMMAL",
        systematic_name="M18009",
        exact_source="hsa04710",
        comment="Sleep disturbance is a core part of ME/CFS",
    )
]
