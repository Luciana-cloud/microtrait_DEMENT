# Y‑A+S is the new Y‑A‑S: Updating microbial life‑history tradeoffs with comparative genomics

Trait‑based ecology proposes that microorganisms balance investments in **yield (Y)**, **resource acquisition (A)**, and **stress tolerance (S)**. These tradeoffs shape microbial community assembly and influence soil biogeochemical processes. This repository contains the data‑processing and analysis workflow for the manuscript:

**_“Updating microbial life history tradeoffs with comparative genomics.”_**

---

## Overview

We analyzed **~30,000 metagenome‑assembled genomes (MAGs)** and **~5,000 soil isolate genomes** to test whether microbial life‑history strategies are encoded at the genome scale. Functional traits were quantified as **gene counts per trait**, and traits were grouped into emergent functional categories aligned with the **YAS framework**:

- **Y (Yield)** – efficiency of converting resources into biomass  
- **A (Acquisition)** – investment in resource uptake and processing  
- **S (Stress tolerance)** – mechanisms enabling survival under environmental stress  

The YAS framework predicts a **three‑way tradeoff**, where investment in one strategy reduces investment in the others.

### Key findings

- No evidence was found for the expected Y–A–S tradeoffs when traits were measured as gene counts.  
- **Genome size** showed a strong positive correlation with **A** and **S** traits.  
- Genome size was **negatively correlated** with **Y**, measured as carbon use efficiency.  
- Larger genomes appear to support broader functional capabilities but incur higher maintenance costs, reducing yield.  
- Genome size may therefore serve as a **master trait** for representing microbial life‑history strategies in trait‑based models.

These results challenge the classical YAS tradeoff model and suggest a revised interpretation: **Y‑A+S**, where acquisition and stress tolerance scale together with genome size, while yield declines.

---

## Repository structure

This repository contains the scripts used to prepare data, generate figures, and reproduce the analyses presented in the manuscript.

1. Final_data_preparation_file.R     # Data preparation and trait aggregation
2. Final_Plots.R                     # Main-text and supplementary figures
3. README.md                         # Project documentation


---

## Workflow summary

### 1. Data preparation  
`Final_data_preparation_file.R` performs:

- Import and cleaning of MAG and isolate metadata  
- Calculation of gene counts per functional trait  
- Aggregation of traits into emergent functional groups  
- Integration of genome size and carbon use efficiency metrics  

### 2. Plot generation  
`Final_Plots.R` produces:

- Genome size vs. trait investment relationships  
- Trait–trait correlation matrices  
- Emergent functional group visualizations  
- All figures used in the main manuscript and supplementary information  

---

## Requirements

The analysis was performed in **R**.  
A typical environment includes:

- `tidyverse`  
- `data.table`  
- `ggplot2`  
- `patchwork`  
- `readr`  
- `dplyr`  

---

## Reproducibility

The scripts are modular and transparent, enabling reuse for:

- Trait‑based ecological analyses  
- Comparative genomics workflows  
- Genome‑scale trait aggregation  
- Model‑ready trait dataset generation  

---

## Contact

For questions or collaboration inquiries, please contact:  
**Luciana Chávez Rodriguez**  
Wageningen University & Research



