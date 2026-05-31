[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.20471842.svg)](https://doi.org/10.5281/zenodo.20471842)

# IOP-dependent and IOP-independent glaucoma component MR code release

This repository provides the public reproducibility code scaffold for a genetic epidemiology study of systolic blood pressure and divergent genetic liability across IOP-dependent and IOP-independent glaucoma components.

## Evidence boundary

The analysis is **hypothesis-generating and non-confirmatory**. The code and manuscript should not be interpreted as establishing a confirmed causal effect of systolic blood pressure on glaucoma, an established IOP-mediated pathway, or patient-level blood-pressure management guidance.

## Repository contents

- `scripts/`: minimal public reproducibility scripts organized by analysis stage.
- `metadata/run_order.tsv`: suggested stage order and example commands.
- `metadata/included_code_manifest.tsv`: included scripts with SHA256 checksums.
- `metadata/excluded_code_manifest.tsv`: files excluded from the minimal release.
- `metadata/release_safety_audit.tsv`: local-path and secret-like pattern audit.
- `config/config.example.yml`: editable example configuration.
- `docs/data_input_requirements.md`: input data requirements and restrictions.

## Data availability

Raw GWAS summary statistics and LD reference panels are not redistributed in this repository. Users must obtain source data from the original providers and comply with their data-use terms.

Generated: 2026-05-24T15:34:15
