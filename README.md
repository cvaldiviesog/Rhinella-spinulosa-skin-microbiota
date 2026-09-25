# Bacterial lineages persist despite microbial network reorganization across amphibian development

This repository contains the R scripts used for the analyses and figure generation presented in:

**Valdivieso C, Chávez FP, Méndez MA, Allende ML.**  
*Bacterial lineages persist despite microbial network reorganization across amphibian development.*

## Repository structure

```text
Rhinella-spinulosa-skin-microbiota/
├── README.md
├── Rhinella_spinulosa_microbiota_pipeline_v2.R
├── Rhinella_spinulosa_microbiota_figures_v2.R

## Scripts

### Rhinella_spinulosa_microbiota_pipeline_v2.R

Main analytical and statistical workflow, including data preprocessing, host-filtering analyses, core microbiota characterization, differential-abundance analyses, microbial association networks, network topology, module and hub analyses, ontogenetic network comparisons, persistence analyses, and sensitivity analyses.

Network inference is performed using SPIEC-EASI, with network analyses and topological characterization incorporating Cytoscape and RCy3.

### Rhinella_spinulosa_microbiota_figures_v2.R

Script for generating the figures presented in the manuscript from the analytical results produced by the main pipeline.
