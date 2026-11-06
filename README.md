# Epilepsy low-VAF SNV analysis

Source data and core analysis code accompanying the epilepsy study.

The study analysed 83 paired epilepsy lesion–peripheral normal
whole-exome sequencing datasets.

## Contents

- Source_data/data/: source data for Figures 2–5 and Supplementary
  Figures 1–5, and a list of 50,000 attempted BAMSurgeon SNV spike-ins,
  including successful-implantation indicators.
- Study_code/sarek/: variant-calling commands and configurations.
- Study_code/analysis/: candidate filtering, fragment-level BAM
  assessment, paired statistics and IGV sampling.
- Study_code/benchmark/: spike-in site screening, implantation,
  high-quality read counting and performance evaluation.

## Requirements

The code uses Python, R, nf-core/sarek and BAMSurgeon.
Python analyses require NumPy and pysam. Resource paths and input
files must be configured before execution.
