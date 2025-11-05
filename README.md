# RiceDBreeder-v1

Dockerized version of **RiceDBreeder** web applications.

---

## 🚀 Quick Start

### 1️⃣ Build & Run

From the project root, run the following command:

```bash
docker compose up -d
```

Then open your browser and access:

```
http://localhost:11028/
```

> ⚙️ If a port conflict occurs, you can edit the port number in `docker-compose.yml` before running.

---

## 🌾 Overview

**RiceDBreeder** integrates three layers of Korean rice data — **pedigree**, **phenotype**, and **genotype** — to enable cross-layer reasoning for digital breeding decisions.

Traditional phenotype-only selection faces challenges in improving multiple traits under climate change, pest pressure, and shifting market demands.
In Korea, phenotypic records, genomic (NGS) data, and pedigree information are separately maintained, limiting their combined use.
RiceDBreeder overcomes these barriers by offering a unified, interactive web interface.

---

## 🏠 Home Interface

The home page provides two primary entry points:

* **Phenotype Search** — for filtering varieties based on phenotypic traits.
* **Rice Variety Search** — for exploring a specific variety’s pedigree, phenotype, and genotype layers.

![Home Page](images/home.png)
*Figure 1. Home interface showing entry points for Phenotype and Variety Search.*

---

## 🧩 Application Modules

### 1️⃣ Phenotype Search

Select **resource type** (*Breeding Line, Cultivar, Germplasm, Landrace, Weedy Type*),
then choose traits such as `planting date`, `grain width`, etc.
The interface provides **distribution plots** and **filter controls** to view and narrow down results.
Filtered varieties are listed in a table with direct links to their integrated variety cards.

![Phenotype Search](images/phenotype_search.png)
*Figure 2. Phenotype Search interface displaying filter panels and trait distribution plots.*

---

### 2️⃣ Rice Variety Search

Serves as the main entry for **Pedigree–Phenotype–Genotype integration**.
Selecting a base variety loads:

* **Pedigree view (left)** — ±2 generations (ancestors and descendants).
* **Analysis panel (right)** — toggles between *Phenotype Analysis* and *GWAS Analysis* tabs.

![Variety Search - Phenotype](images/variety_search_pheno.png)
*Figure 3. Rice Variety Search — Phenotype analysis showing trait distributions with the base variety highlighted.*

![Variety Search - Genotype](images/variety_search_geno.png)
*Figure 4. Rice Variety Search — Genotype analysis showing subtrait-based variant counts and GWAS scatter plots.*

---

## 🌳 Pedigree Visualization

* **Scope:** Two parental generations up and two descendant generations down.
* **Node colors:**

  * 🟦 Genotype only
  * 🟩 Phenotype only
  * 🟧 Both genotype & phenotype
  * ⬜ No linked data
* **Shapes:**

  * Round = found in pedigree
  * Rectangle = not connected to pedigree

### Pedigree Features

* **Expand / Remove:** Left-click to expand or collapse one generation.
* **Reset:** Return to the initial state.
* **Add Additional Varieties:**

  * Pedigree-based search
  * Phenotype-based search (same filters as Phenotype Search)
  * Genotype-based search

---

## 🧬 Genotype Analysis & Filtering Options

RiceDBreeder provides two integrated filtering systems — **Filter Options** and **Group Options** — both of which support trait/subtrait selection.

### 🧠 1. Filter Options (with Subtrait Filtering)

| Option                       | Description                                                                                   |
| ---------------------------- | --------------------------------------------------------------------------------------------- |
| **Trait/Subtrait Selection** | Select domains (e.g., yield, stress, biochemical). Determines which variants appear in plots. |
| **P-value**                  | Adjust GWAS significance threshold (red line in scatter plot).                                |
| **MAF Filter**               | Minimum minor allele frequency (default ≥ 0.05).                                              |
| **SNP Presence Threshold**   | Filter variants appearing in ≥ N samples.                                                     |
| **Unique Mode**              | Show sample-unique variants (2–5 samples).                                                    |

---

### 🧩 2. Group Options (with Trait-Based Selection)

| Option                            | Description                                                               |
| --------------------------------- | ------------------------------------------------------------------------- |
| **Group Selection (Trait-based)** | Create up to 5 comparison groups based on selected traits or subtraits.   |
| **Variant ID Filter**             | Compare results between different presence thresholds (e.g., N=4 vs N=6). |
| **Sample Grouping**               | Identify variants unique to each selected sample group.                   |

Trait-driven grouping enables comparison of allele distribution across categories (e.g., yield vs stress).

---

## 🔄 Pedigree–Analysis Interactions

* **Phenotype tab:** Clicking nodes highlights phenotype markers.
* **Genotype tab:** Clicking nodes adds pink borders and updates bar/scatter overlays.
* Interactions are **bidirectional**, ensuring that all panels remain synchronized.

---

## ⚙️ Features Summary

* Expandable pedigree (+1 generation per click)
* Multi-variety integration
* Integrated trait-based filtering across all views
* Bidirectional linkage between graph and plots
* Reset, remove, and memory-preserved filters
* GWAS–VCF–Pedigree synchronization
