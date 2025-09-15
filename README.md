## St. Lawrence Connectivity

This repository contains analyses of ecological connectivity prioritization in the St. Lawrence Lowlands ecoregion in Quebec, Canada.
Results are summarized in two technical reports and in a manuscript currently under review.

Technical report — Present & future results (French)
https://www.environnement.gouv.qc.ca/biodiversite/cadre-ecologique/modelisation-connectivite-basses-terres-saint-laurent.pdf

Technical report — Present-only results (French)
https://www.environnement.gouv.qc.ca/biodiversite/cadre-ecologique/priorisation-connectivite-basses-terres-saint-laurent.pdf

Manuscript in review: Linking multispecies connectivity models and long-term scenarios to guide conservation — Rayfield, Boulanger, Larocque, Lucet, Teixeira-Martins, and Gonzalez (in review).

## Repository Structure

```text
st-lawrence-connectivity/
├─ present-and-future prioritization/      # Mixed scripts: present + future scenarios across regions
│  ├─ data/                                # Inputs (rasters, species lists, masks, params)
│  ├─ scripts/                             # Analysis code reused across regions (R, shell)
│  ├─ outputs/                             # Intermediate/final results (tables, rasters)
│  └─ figures/                             # Publication-ready plots/maps
│
├─ present-only prioritization/            # Present-day runs by CERQ Level-2 “natural regions”
│  ├─ b01b02/                              # B01+B02 combined (Plaine du haut & moyen Saint-Laurent)
│  │   ├─ data/
│  │   ├─ scripts/
│  │   ├─ outputs/
│  │   └─ figures/
│  └─ b03/                                 # B03 (Plaine d’Ottawa)
│      ├─ data/
│      ├─ scripts/
│      ├─ outputs/
│      └─ figures/
│
├─ .gitignore
└─ README.md
