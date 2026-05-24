# Public code release notes

This repository contains a minimal public reproducibility code release for a summary-statistics genetic epidemiology analysis of systolic blood pressure and IOP-dependent versus IOP-independent glaucoma genetic components.

## Claim boundary

The code release supports a hypothesis-generating and non-confirmatory analysis framework. It should not be interpreted as establishing a confirmed causal effect of systolic blood pressure on glaucoma, a confirmed IOP-mediated pathway, or patient-level blood-pressure management guidance.

## Data boundary

Raw GWAS summary statistics and restricted third-party data are not redistributed in this repository. Users should obtain source GWAS files from the original providers and configure local paths using `config/config.example.yml`.

## Recommended use

Review `metadata/run_order.tsv`, `metadata/data_availability.tsv`, and `docs/data_input_requirements.md` before attempting to rerun the workflow.
