
# CMV in lung transplant recipients: analysis code

R code used for the statistical analysis and figure generation of the doctoral
thesis *Aspectos virológicos de la infección por citomegalovirus en receptores
de trasplante pulmonar* (Universidad Autónoma de Madrid, 2026).

Each code section is labelled with the corresponding figure or table number in
the thesis and with the statistical test applied, so that every reported result
can be traced back to the code that produced it.

## Contents

- `study1-genotyping/` — Retrospective cohort of lung transplant recipients
  (2009–2021). Genotyping of CMV glycoproteins gB (UL55) and gH (UL75) by
  genotype-specific real-time PCR. Analyses at sample, episode and patient level.
- `study2-ngs-resistance/` — Characterisation of antiviral resistance variants
  in UL97, UL54 and UL56 by next-generation sequencing (DeepChek CMV v2.0),
  classified according to the CHARMD framework.

## Data availability

**This repository contains code only.** No patient-level data are included.

The scripts expect data files that are not distributed here for confidentiality
reasons. The expected structure of each input file is described in
`data/README.md`. Access to the underlying data may be requested through the
corresponding author, subject to institutional and ethical approval.

## Requirements

R (≥ 4.3) and the following packages: `tidyverse`, `survival`, `flextable`,
`officer`, `ggplot2`, `patchwork`, `readxl`. Exact versions are listed in
`sessionInfo.txt`.

## How to run

Scripts are organised in foldable RStudio sections marked as
`## ---- [SECTION] ----`. Each section is self-contained and can be executed
independently once the data objects have been loaded by the setup section at
the top of each file.

## Citation

Reyes Ruiz CA. *Aspectos virológicos de la infección por citomegalovirus en
receptores de trasplante pulmonar* [doctoral thesis]. Madrid: Universidad
Autónoma de Madrid; 2026.

## Author

Carmen Alhena Reyes Ruiz — ORCID [0000-0003-0778-4971](https://orcid.org/0000-0003-0778-4971)

Thesis supervised by Dra. María Dolores Folgueira López, Servicio de
Microbiología, Hospital Universitario 12 de Octubre, Madrid.

## License

Released under the MIT License. See `LICENSE`.
