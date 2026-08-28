from typing import get_args

from mecfs_bio.assets.reference_data.polyfun.annotations.baseline_lf_annotation_names import (
    BASELINE_LF_ANNOTATION_NAMES,
)
from mecfs_bio.constants.polyfun_annotation_families import (
    FAMILY_SHORT_LABELS,
    AnnotationFamily,
    family_for_annotation,
)


def test_names_constant_is_187_unique():
    assert len(BASELINE_LF_ANNOTATION_NAMES) == 187
    assert len(set(BASELINE_LF_ANNOTATION_NAMES)) == 187


def test_every_annotation_maps_to_a_valid_family():
    valid = set(get_args(AnnotationFamily))
    for name in BASELINE_LF_ANNOTATION_NAMES:
        assert family_for_annotation(name) in valid


def test_representative_family_assignments():
    cases = {
        "non_synonymous_lowfreq": "non_synonymous",
        "Coding_UCSC_common": "coding",
        "synonymous_common": "coding",
        "UTR_3_UCSC_common": "coding",
        "Conserved_Primate_phastCons46way_common": "conserved",
        "GERP.RSsup4_common": "conserved",
        "GERP.NS_common": "conserved",
        "Promoter_UCSC_common": "promoter_or_enhancer",
        "TSS_Hoffman_common": "promoter_or_enhancer",
        "SuperEnhancer_Hnisz_common": "promoter_or_enhancer",
        "BivFlnk.flanking.500_common": "promoter_or_enhancer",
        "H3K27ac_Hnisz_common": "histone_marks",
        "H3K4me1_peaks_Trynka_common": "histone_marks",
        "Repressed_Hoffman_common": "repressed",
        "DHS_Trynka_common": "open_chromatin",
        "DHS_peaks_Trynka_common": "open_chromatin",
        "FetalDHS_Trynka_common": "open_chromatin",
        "DGF_ENCODE_common": "open_chromatin",
        "TFBS_ENCODE_common": "other",
        "CTCF_Hoffman_common": "other",
        "Transcr_Hoffman_common": "other",
        "Intron_UCSC_common": "other",
        "MAFbin_frequent_3": "maf_bins",
        "CpG_Content_50kb_common": "ld_related_continuous",
        "Recomb_Rate_10kb_common": "ld_related_continuous",
        "Backgrd_Selection_Stat_common": "ld_related_continuous",
        "MAF_Adj_ASMC_common": "ld_related_continuous",
        "GTEx_eQTL_MaxCPP_common": "molecular_qtl",
        "BLUEPRINT_H3K27acQTL_MaxCPP_common": "molecular_qtl",
    }
    for name, family in cases.items():
        assert family_for_annotation(name) == family


def test_short_labels_cover_every_family():
    assert set(FAMILY_SHORT_LABELS) == set(get_args(AnnotationFamily))


def test_open_chromatin_membership_is_exactly_dhs_and_dgf():
    open_chrom = {
        n
        for n in BASELINE_LF_ANNOTATION_NAMES
        if family_for_annotation(n) == "open_chromatin"
    }
    # every open_chromatin member is a DHS/FetalDHS/DGF accessibility annotation
    assert all(("DHS" in n) or ("DGF" in n) for n in open_chrom)
    # and every DHS/DGF accessibility annotation is captured
    assert all(
        family_for_annotation(n) == "open_chromatin"
        for n in BASELINE_LF_ANNOTATION_NAMES
        if ("DHS" in n) or ("DGF" in n)
    )
