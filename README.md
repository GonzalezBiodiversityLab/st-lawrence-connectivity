## St. Lawrence Connectivity

This repository contains analyses of ecological connectivity prioritization in the St. Lawrence Lowlands ecoregion in Quebec, Canada.
Results are summarized in two technical reports and in a manuscript currently under review.

Technical report — Present & future results (French)
https://www.environnement.gouv.qc.ca/biodiversite/cadre-ecologique/modelisation-connectivite-basses-terres-saint-laurent.pdf

Technical report — Present-only results (French)
https://www.environnement.gouv.qc.ca/biodiversite/cadre-ecologique/priorisation-connectivite-basses-terres-saint-laurent.pdf

Manuscript under review: Linking multispecies connectivity models and long-term scenarios to guide conservation — Rayfield, Boulanger, Larocque, Lucet, Teixeira-Martins, and Gonzalez (in review).

## Repository Structure

```text
st-lawrence-connectivity/
├─ present-and-future prioritization/      # Scenario-based analyses (present + future) across CERQ natural regions
│  ├─ data/                                # Inputs (rasters, species lists, masks, params)
│  ├─ scripts/                             # Analysis code shared across regions (R/shell)
│  └─ libraries/                           # Parameterized ST-Sim library of landscape change
└─ present-only prioritization/            # Present-day analyses by CERQ Level-2 natural regions (régions naturelles: B01, B02, B03)
   ├─ b01b02/                              # B01+B02 combined (Plaine du haut & moyen Saint-Laurent)
   │  └─ scripts/                          # Subregion-specific code for the combined B01+B02 runs
   └─ b03/                                 # B03 (Plaine d'Ottawa)
      ├─ data/                             # Inputs specific to B03 adapted from B02
      └─ scripts/                          # Subregion-specific code for B03

