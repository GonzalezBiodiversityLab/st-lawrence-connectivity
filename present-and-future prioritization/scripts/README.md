# Scripts - present-and-future prioritization (archival)

**Purpose.** This folder contains an archival snapshot of analysis scripts used by four analysts for the “present-and-future prioritization” work. Scripts are provided **as-is** in the spirit of open science and transparency.

**Important notes**
- These scripts are **not a unified, runnable pipeline**. Different files assume different directory structures and machine-specific paths.
- We have **not** modified scripts post-analysis to avoid introducing discrepancies with the published results.
- Re-running any script will likely require editing file paths and having access to the original input data, which are not all publicly hosted.

**What’s here**
- Individual R scripts and utilities used during data preparation, modeling, and post-processing.
- See the **script index** below for a one-line description of each file.

**Compute environment (high level)**
- Primary language: **R** (with common spatial/analysis packages).
- Other tooling used in the project: **Circuitscape**, **GDAL/QGIS**, **SyncroSim**.
- OS: mixed (analyst-dependent).

**Known limitations**
- Analyst-specific absolute paths.
- Heterogeneous style and assumptions across scripts.
- Some scripts reference inputs that are large or not publicly hosted.

---

## Script index (archival)

_Scripts are provided as-is from multiple analysts; paths and assumptions differ. The index is ordered by the leading step numbers (e.g., 09-01, 09-02). These numbers are for navigation only; scripts were authored independently._

| Step | File | Brief purpose |
|---:|---|---|
| 00 | `00_constants.R` | Define shared constants, options, and reusable paths/parameters for downstream scripts 01 - 07. |
| 01 | `01_crop-strata.R` | Prepare analysis strata (e.g., natural region strata, protected areas) used in landscape change model. |
| 02 | `02_mrc-area.R` | Compute MRC-level area summaries per stratum. |
| 03 | `03_state-class-forest-composition-and-age.R` | Assemble state-class inputs for forest composition and age initial conditions. |
| 04 | `04_transition-multipliers.R` | Build transition multiplier layers for climate scenarios. |
| 05 | `05_protected-areas-spatial-multipliers.R` | Create spatial multipliers reflecting protected areas. |
| 06 | `06_model-post-processing-crop-rasters.R` | Post-process model landscape rasters to crop to study area. |
| 07 | `07_postprocess-prob-conversion.R` | Derive probability-of-conversion rasters and summaries. |
| 08 | `08-01_postprocess-habitat-connectivity-percent-change.R` | Calculate percent change in habitat/connectivity between baseline and scenarios. |
| 08 | `08-02_plot-habitat-connectivity-percent-change.R` | Plot percent-change connectivity results for reporting. |
| 09 | `09-01_connectivity-DataExtractionFormatting.R` | Organize data from connectivity analyses (Circuitscape). |
| 09 | `09-02_connectivity-MappingCurrentDensity.R` | Map current density results from connectivity analyses. |
| 09 | `09-03_connectivity-MappingEcoregionConnectivity.R` | Map/aggregate connectivity metrics by ecoregion/natural region. |
| 09 | `09-04_connectivity-GraphicalSummaries.R` | Produce graphical summaries (figures/plots) of connectivity results. |
| 09 | `09-05_connectivity-CompositeFigures.R` | Assemble composite connectivity figures for the report/manuscript. |
| 10 | `10_zonation-preprocess.R` | Prepare Zonation prioritization inputs. |
| 11 | `zonation_run/` | Shell scripts (.sh) to prepare and run Zonation analyses on the Béluga HPC cluster (Digital Research Alliance of Canada). |


> This index is descriptive only; it does not guarantee re-runnability.
