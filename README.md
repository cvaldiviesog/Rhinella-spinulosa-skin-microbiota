# Host-associated filtering and microbial network organization in the skin microbiota of an Andean amphibian

This repository contains the R scripts used for the analyses and figure generation presented in:

**Valdivieso C, Chávez FP, Méndez MA, Allende ML.**  
*Host-associated filtering and microbial network organization in the skin microbiota of an Andean amphibian.*

## Repository structure

```text
Rhinella-spinulosa-skin-microbiota/
├── README.md
├── Rhinella_spinulosa_microbiota_pipeline_v2.R
└── Rhinella_spinulosa_microbiota_figures_v2.R
```

## Scripts

### `Rhinella_spinulosa_microbiota_pipeline_v2.R`

Main analytical and statistical workflow, including data preprocessing, host-associated filtering analyses, core microbiota characterization, differential-abundance analyses, microbial association network inference, network topology, module and hub analyses, sensitivity analyses, and sequence-level correspondence across sampling groups.

Network inference is performed using SPIEC-EASI. Network visualization and topological characterization incorporate Cytoscape and RCy3.

### `Rhinella_spinulosa_microbiota_figures_v2.R`

Script for generating the figures presented in the manuscript from the analytical results produced by the main pipeline.
