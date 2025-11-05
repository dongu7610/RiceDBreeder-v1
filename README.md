# RiceDBreeder-v1

Dockerized version of **RiceDBreeder** web applications.

---

## 🚀 Quick Start

### 1️⃣ Build & Run

From the project root, run the following command:

```bash
docker compose up -d
```



> ⚙️ If a port conflict occurs, you can edit the port number in `docker-compose.yml` before running.

---

## 🌾 Overview

**RiceDBreeder** integrates three layers of Korean rice data — **pedigree**, **phenotype**, and **genotype** — to enable cross-layer reasoning for digital breeding decisions.

Traditional phenotype-only selection faces challenges in improving multiple traits under climate change, pest pressure, and shifting market demands.
In Korea, phenotypic records, genomic (NGS) data, and pedigree information are separately maintained, limiting their combined use.
RiceDBreeder overcomes these barriers by offering a unified, interactive web interface.

---

## 🧩 Application Modules

### 1️⃣ Phenotype Search

* Select **resource type**: *Breeding Line, Cultivar, Germplasm, Landrace, Weedy Type*.
* Choose traits such as `planting date`, `grain width`, etc.
* Displays **distribution plots** showing full variety-level variation with filter controls.
* Filtered results appear as an interactive **variety table**.
* Clicking a variety ID (or using the search bar in **Rice Variety Search**) opens the integrated analysis page.

### 2️⃣ Rice Variety Search

* Serves as the main entry for **Pedigree–Phenotype–Genotype integration**.
* Selecting a base variety loads:

  * **Pedigree view (left)** — ±2 generations (ancestors and descendants).
  * **Analysis panel (right)** — toggles between *Phenotype* and *Genotype* tabs.

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

If multiple pedigrees overlap, they merge into a single unified graph; if disconnected, they remain separate.
Expanded nodes (orange/green/blue) can be clicked to sync phenotype/genotype panels — highlighted in **pink** on selection.

---

## 📊 Phenotype Analysis Tab

* Displays global trait distributions (e.g., planting date, grain width).
* Highlights the searched variety within each plot.
* Clicking a trait sends its values to the pedigree, coloring nodes accordingly.

---

## 🧬 Genotype Analysis Tab

Integrates **VCF genotype data** and **GWAS marker–trait associations**.

* Only variants where the searched variety carries the **GWAS minor allele** are shown.
* Each variant is linked to a trait and visualized in:

  * **Trait-count bar plot** (colored by subtrait group)
  * **Scatter plot (Chr1–Chr12)**
  * **Variant data table**

Clicking bar or scatter elements highlights corresponding pedigree nodes.

---

## 🧮 Filtering Options

RiceDBreeder provides three main filtering layers:

### 🧠 1. Filter Options

| Option                     | Description                                                 |
| -------------------------- | ----------------------------------------------------------- |
| **P-value**                | Adjust GWAS significance cutoff (red line in scatter plot). |
| **MAF Filter**             | Minimum minor allele frequency threshold (default ≥ 0.05).  |
| **SNP Presence Threshold** | Filter variants appearing in ≥ N samples.                   |
| **Unique Mode**            | Show sample-unique variants (for 2–5 samples).              |

### 🧩 2. Group Options

| Option                | Description                                                               |
| --------------------- | ------------------------------------------------------------------------- |
| **Group Selection**   | Create up to 5 comparison groups (by trait, variant, or sample).          |
| **Variant ID Filter** | Compare results between different presence thresholds (e.g., N=4 vs N=6). |
| **Sample Grouping**   | Identify variants unique to each group of selected samples.               |

### 🌿 3. Trait Options

| Option                       | Description                                              |
| ---------------------------- | -------------------------------------------------------- |
| **Trait/Subtrait Selection** | Choose trait domains (e.g., yield, stress, biochemical). |
| **Trait-based Filtering**    | Restrict analysis to specific subtrait-linked variants.  |

---

## 🔄 Pedigree–Analysis Interactions

* **Phenotype tab:** Clicking nodes highlights phenotype markers.
* **Genotype tab:** Clicking nodes adds pink borders and updates bar/scatter overlays.
* Interactions are **bidirectional** — selections in one view update all others.

---

## ⚙️ Features Summary

* Expandable pedigree (+1 generation per click)
* Multi-variety integration
* Interactive filtering (trait, variant, group)
* Bidirectional linkage between graph and plots
* Reset, remove, and filter memory functions

---



