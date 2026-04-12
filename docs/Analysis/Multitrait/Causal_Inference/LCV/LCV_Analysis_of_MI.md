#  Analysis of MI

To investigate the reliability of the Latent Causal Variable (LCV)[@o2018distinguishing] causal inference technique, I used LCV to estimate the Genetic Causality Proportion (GCP) of several well-known risk factors on the risk of myocardial infarction (MI). I also included educational attainment as a negative control.  The results are below.


{{ include_file("docs/_figs/mi_analysis_lcv_table_downstream_MI.mdx") }}


- The strong effect of LDL and triglycerides on risk of MI is consistent with their extensively documented status as causal risk factors[@steinberg2011cholesterol].
- The high GCP for C-reactive Protein (CRP) is interesting.  A cursory search suggests that while CRP is a very strong biomarker for inflammation, which is in turn a strong risk for cardiovascular disease, there is controversy about whether CRP itself actually plays a causal role.  Thus The high GCP score may indicate a limitation of LCV as a causal inference technique.
- The low GCP for educational attainment indicates it has successfully served as a negative control.

