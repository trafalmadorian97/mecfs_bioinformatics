# DecodeME preprint and data analysis plan

Primary sources for the question of whether DecodeME build-38 summary statistics
are reference-aligned. The PDFs are gitignored (*.pdf); the extracted text
(*.txt, via pdftotext -layout in the pdf-env) is tracked.

Sources (medRxiv itself 403s to scripted clients; the Edinburgh mirror serves
the PDFs):
- Preprint: medRxiv 10.1101/2025.08.06.25333109v1
  https://institute-genetics-cancer.ed.ac.uk/sites/default/files/2026-05/2025-08-03%20DecodeME%20Preprint.pdf
- Supplementary Methods: provided by the user (decode_me_supplement.pdf); the
  medRxiv supplement for 10.1101/2025.08.06.25333109.
- Data Analysis Plan (2024):
  https://institute-genetics-cancer.ed.ac.uk/sites/default/files/2025-06/DecodeME%20data%20analysis%20plan.pdf
- Summary statistics (GWAS-1): OSF project rgqs3.

Re-fetch and extract:
  UA="Mozilla/5.0 (X11; Linux x86_64) AppleWebKit/537.36"
  curl -sL -A "$UA" -o decodeme_preprint.pdf "<preprint url above>"
  curl -sL -A "$UA" -o decodeme_data_analysis_plan.pdf "<DAP url above>"
  pixi r -e pdf-env pdftotext -layout decodeme_preprint.pdf decodeme_preprint.txt
