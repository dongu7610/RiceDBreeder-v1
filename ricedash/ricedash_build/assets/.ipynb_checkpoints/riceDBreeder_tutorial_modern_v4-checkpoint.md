<!-- Final: 3-section layout (background / preview / using) -->

<h2 id="top"></h2>

# RiceDBreeder · Tutorial

*A compact guide to why RiceDBreeder exists, what it integrates, and how to use it.*

---

<h2 id="background">1. Background & Data Integration</h2>

Traditional, phenotype-only selection struggles with **multi-trait improvement** under climate change, emerging pests/diseases, and shifting market demands.
In Korean rice resources, integration is further hindered because phenotype and genomics are stored separately, **NGS data is fragmented**, and pedigrees remain in unstructured documents — making cross-layer reasoning difficult.

**What we integrate**

* **Phenotype (Genebank):** 8,831 varieties, 62 standardized traits.
* **Pedigree (RDA):** 843 curated varieties converted from HWP/PDF lineage records.
* **Resequencing (NABIC VCF):** SNP data indexed per variety for fast retrieval.
* **GWAS Atlas:** 163,479 marker–trait associations across 461 traits (475 resequenced varieties).

**Why it matters**

* Run **multi-trait filtering** on unified phenotype data.
* Inspect **trait-linked SNPs** per variety and compare among groups.
* Use **pedigree context** to trace inheritance and relatedness.
* Move from fragmented spreadsheets to a **coherent digital breeding workflow**.

**Objectives**

1. **Integrate** phenotype, pedigree, GWAS, and resequencing data.
2. **Deliver** an interactive Plotly Dash application for filtering and exploration.
3. **Enable** rapid cross-linked navigation for breeding insights.

---

<h2 id="preview">2. Overview Preview</h2>

RiceDBreeder connects three layers — phenotype, pedigree, and genotype — into a single, interactive system.

![Web overview](/assets/web_preview.png)

---

<h2 id="using">3. Using the App</h2>

RiceDBreeder provides two primary pages:

* **Phenotype Search** — Filter and select varieties by traits.
* **Integrated Visual Dashboard** — Explore pedigree, phenotype, and genotype information interactively.

---

<h3 id="home">Home Interface</h3>

From the home page, users can start from either **Phenotype Search** or **Rice Variety Search**.

* **Phenotype Search:** Explore accessions by phenotype filters.
* **Rice Variety Search:** Jump directly to a specific variety card containing pedigree and genomic information.

![Home interface](/assets/home.png)

---

<h3 id="phenotype-search">Phenotype Search</h3>

Select a **resource type** such as “Breeding Line”, “Cultivar”, “Landrace”, or “Weedy Type”.
Choose phenotype traits (e.g., `planting date`, `grain width`) and define filter ranges.
The app visualizes trait distributions and displays a **filtered variety list**.
Each listed variety includes an **IT number** link to access its detailed card.

![Phenotype Search page](/assets/page1.png)

---

<h3 id="visual-dashboard">Integrated Visual Dashboard</h3>

When a variety is opened from **Rice Variety Search** or **Phenotype Search**,
it becomes the **base variety** for analysis.
The dashboard combines:

* **Pedigree Visualization (left)** — An expandable ±2-generation network centered on the base variety.
* **Analysis Dashboard (right)** — Switchable **Phenotype** and **GWAS/SNP** tabs.

**Pedigree Visualization**

* Node color indicates data type:

  * 🟦 Genotype only
  * 🟩 Phenotype only
  * 🟧 Both
  * ⬜ No data
* Node shape:

  * **Circle:** in pedigree
  * **Rectangle:** no pedigree record

---

<h2 id="controls">Pedigree Interaction Controls</h2>

**Core controls:**

* **Expand:** Add one more generation of parents or offspring.
* **Remove:** Delete selected nodes except the base.
* **Reset:** Return to the initial layout including additions.
* **Wide:** Expand the canvas width for complex pedigrees.

---

<h2 id="add-option">Add Option Overview</h2>

Add varieties to the pedigree using three methods:

1. **Pedigree-based Add:** Search by name and expand its ±2-generation structure.
2. **Phenotype-based Add:** Select varieties meeting Page 1–style filters.

   * Linked (circle) if pedigree exists
   * Standalone (rectangle) if missing
3. **Genotype-based Add:** Add varieties with genotype data for joint analysis.

After selection, press **Apply** to merge new pedigrees with existing ones.
When multiple pedigrees overlap, they merge automatically.

![Highlighted pedigrees](/assets/pedigree1_pedigree2.png)

---

<h2 id="phenotype">Phenotype Tab</h2>

Links the **pedigree (left)** and **phenotype view (right)**.
Selecting nodes with phenotype data (🟩 *Phenotype only*, 🟧 *Both*)
shows their position in the full trait distribution.

* Click nodes to visualize selected traits in histogram/barplot form.
* The same color gradient reflects back to pedigree nodes for consistency.
* Use this to compare relatives’ phenotypic positions efficiently.

---

<h2 id="gwas">GWAS Tab</h2>

For varieties with genotype data (🟦 *Genotype only*, 🟧 *Both*),
the GWAS tab enables filtering and visualization of associated SNPs.

### Filter Options

Adjust filtering by **Subtrait**, **P-value cutoff**, **MAF threshold**, **SNP presence**, or **Unique Mode**.

![Filter options](/assets/filter_option.png)

### Pedigree Integration Examples

Static figures illustrate how selections and filters affect the pedigree view:

* **Variant selection highlight:**
  ![Selected variant](/assets/select_variant.png)

* **Variant reflection to pedigree:**
  ![Variant reflection](/assets/select_variant_result.png)

* **Unique mode active:**
  ![Unique scatter result](/assets/unique_gwas_scatter_plot.png)

* **No-pedigree varieties after Add Option:**
  ![Additional nopedi case](/assets/additional_nopedi.png)

---

<h2 id="top-link"></h2>
