# Genomic data conventions

Rules about the meaning and provenance of biological data, beyond code style.

## Column semantics

- **A column holds exactly what its name says.** When the true value can't be determined, write
  null. Never put a stand-in that happens to give the right answer for today's consumer (e.g. minor
  allele frequency in an effect-allele-frequency column). Check every consumer of the column, and
  make sure they handle nulls.

## Gene and protein identifiers

- **Never join tables on gene symbols.** Symbols go stale and collide across loci and proteins.
  Join on stable ids (UniProt accession, Ensembl gene id, HGNC id, Entrez id). When two tables share
  no id, bridge them through an id-mapping asset such as the HGNC complete set. Take display
  symbols from the id-keyed table, not from the source. A few unmapped rows (null plus a warning)
  are better than a symbol fallback.
- **Get gene and protein information from in-repo assets, not from gget or per-item web APIs.**
  Search `mecfs_bio/assets/reference_data/` first. Examples:
  - `HGNC_UNIPROT_GENES` in `hgnc/hgnc_complete_set.py` maps UniProt accessions to genes.
  - `UNIPROT_LOOKUP` in `uniprot/uniprot_lookup_table.py` holds protein names, function, and
    subcellular location.
  - `GENE_THESAURUS` in `ensembl_biomart/gene_thesaurus.py` holds Ensembl ids, symbols, and
    coordinates; join it on Ensembl id.

  Don't add new gget call sites.

## Allele harmonization

- **`GenomeReferenceHarmonizationTask`
  (`mecfs_bio/build_system/task/genome_reference_harmonization/`) is the harmonizer.** It orients
  each variant so that NEA is the plus-strand reference allele and flips allele-dependent
  statistics. Before adding a harmonization step, check whether this Task already runs upstream.
- **Name the harmonization step you mean.** "Harmonization" covers distinct operations:
  - genome-reference harmonization (the Task above)
  - gwaslab harmonization (gwaslab's `harmonize`)
  - task-reference harmonization (an allele join against a reference table, e.g.
    `HarmonizeGWASWithReferenceViaAlleles`)

  Use the specific name in prose, specs, and docstrings.

## Significance thresholds

- **Protein-side thresholds stay genome-wide significant (p ≤ 5e-8).** This covers pQTL hit tables
  and protein screening. Relaxed thresholds (e.g. 1e-5) are for exploring weak signals on the
  trait side only.
