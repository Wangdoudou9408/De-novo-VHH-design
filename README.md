# NK2R VHH Agonist Design Pipeline

De novo design pipeline for VHH nanobodies targeting the Neurokinin 2 Receptor (NK2R) activation pocket, developed for [PDC2026 Competition](https://lab.cs.tsinghua.edu.cn/scc/competition/PDC2026).

## Overview

This repository contains the complete computational workflow for designing VHH binders that activate NK2R by mimicking Neurokinin A (NKA) binding to the orthosteric pocket.

## Pipeline Steps

1. **RFdiffusion**: Generate VHH backbone scaffolds targeting NK2R hotspot residues
2. **Tier 1 Filtering**: Geometric/interface quality assessment
3. **ProteinMPNN**: Sequence design for generated backbones
4. **Tier 2 Scoring**: Designability assessment
5. **Tier 3 Filtering**: Expression quality metrics (GRAVY, pI, aggregation, disulfide bonds)
6. **Boltz-1/AF3**: Structure prediction and interface validation
7. **GROMACS MD**: Stability verification (optional)

## Installation

```bash
# Clone repository
git clone https://github.com/Wangdoudou9408/nk2r-vhh-design.git
cd nk2r-vhh-design

# Install dependencies (example)
conda create -n vhh_design python=3.12
conda activate vhh_design
pip install biopython numpy scipy pandas matplotlib
# Add other dependencies as needed
