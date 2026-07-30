# deCODE rheumatoid arthritis cohort sample sizes

`decode_cohort_sample_sizes.csv` records the per-cohort case/control counts for the
deCODE rheumatoid arthritis meta-analyses, transcribed from **Table 1** of:

> Saevarsdottir S, Stefansdottir L, Sulem P, et al. "Multiomics analysis of
> rheumatoid arthritis yields sequence variants that have large effects on risk of
> the seropositive subset." *Annals of the Rheumatic Diseases* 2022;81:1085-1095.

The deCODE summary statistics do not carry a sample-size column. Instead each
variant has a `Cohorts` string whose characters (`+`/`-`/`?`, where `?` means the
cohort did not contribute to that variant) report, in a fixed order, which cohorts
contributed. Per the README shipped with the download, that order is:

    Iceland, Norway, Sweden, Denmark, UK, Finland

The `position` column gives the 0-based index of each cohort within the `Cohorts`
string, so the table can be used to turn the per-variant membership string into a
per-variant effective sample size (sum over contributing cohorts of
`4 / (1/n_cases + 1/n_controls)`). See
`mecfs_bio/build_system/task/pipes/effective_n_from_cohort_string_pipe.py`.

Both serostatuses are included so the seronegative analysis can reuse this table.
Row totals reconcile with the paper: seropositive = 18,019 cases / 991,604 controls;
seronegative = 8,515 cases / 1,015,471 controls.
