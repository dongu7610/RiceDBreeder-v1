<!-- Final: 3-section layout (background / preview / using) -->
<h2 id="top"></h2>

# RiceDBreeder · Tutorial

_A compact guide to why RiceDBreeder exists, what it integrates, and how to use it._

---

<h2 id="background">1. Background & Data Integration</h2>


Traditional, phenotype-only selection struggles with **multi-trait** improvement under climate change, emerging pests/diseases, and shifting market demands.  
In Korean rice resources, integration is further hindered because phenotype and genomics are stored separately, **NGS data is fragmented**, and pedigrees remain in unstructured documents—making cross-layer reasoning difficult.

**What we integrate**
- **Phenotype (Genebank):** 8,831 varieties, **62 traits** (standardized).  
- **Pedigree (RDA docs):** **843 varieties** curated from HWP/PDF into tidy parent–child tables.  
- **Resequencing (NABIC VCF):** per-variety SNPs indexed for fast retrieval.  
- **GWAS Atlas:** **163,479** marker–trait associations across **461 traits**; minor-allele-linked SNPs condensed per variety (**475 varieties**).

**Why it matters**
- Run **multi-trait** queries on clean phenotype data (not isolated spreadsheets).  
- Inspect **trait-linked SNPs** per variety and compare across groups.  
- Use **pedigree context** to explore related varieties and trace signals end-to-end.  
- Move from ad-hoc lookups to a **unified, decision-ready** workflow.

**Objectives**
1) **Integrate** phenotype, pedigree, GWAS, resequencing into one navigable layer.  
2) **Deliver** a Plotly Dash app enabling filter → inspect → compare.  
3) **Enable** fast visualization and cross-navigation for breeding decisions.

---

<h2 id="preview">2. Overview Preview</h2>

The figure below summarizes the integrated database layers and the two main pages.

![Web preview](/assets/web_preview.png)

---

<h2 id="using">3. Using the App</h2>

**App structure**  
- `home.py` — entry with links to **Tutorial** and modules.  
- `app7_1.py` — Page 1 (multi-filter for phenotypes).  
- `app7_2.py` — Page 2 (resource detail: basic info, **Pedigree**, **Phenotype**, **GWAS/SNP**).

---
<h2 id="page1">Page 1 · Multi-filter (phenotype)</h2>

- Select traits and **stack filters** (e.g., A → 2,000; A + B → 500).  
- Produce a **result table** and pass selection to Page 2.  

In the **Phenotype Search** interface,  
select resource type: _Breeding Line_, adjust the **planting** and **seeding date** ranges,  
and review key results — **IT number**, **VCF**, **pedigree**, and **phenotype values**.

![Page 1 – Multi-filter](/assets/page1.png)

---
<h2 id="page2">Page 2 · Integrated Visual Dashboard</h2>

Arrive here from **Home → Page 2** or **Page 1 → Page 2**.  
The variety you used to enter becomes the **base variety**.

- **Pedigree**: Draws a ±2-generation graph centered on the base variety (expandable on demand).  
  Use **Add Option** to bring in additional context (e.g., candidates not in the current pedigree) for linked analysis.

- **Phenotype**: For each trait, shows the full distribution and clearly marks the **base variety’s position** (helpful for quick percentile/context checks).

- **GWAS/SNP**: Displays **GWAS-related SNP positions across chr1–chr12** for the base variety, with overlay/filters to inspect relevant signals.
---
<h2 id="pedigree">Pedigree Overview</h2>

Pedigree graph visualizes family structure of the selected variety.  
Node colors indicate data availability status:

- **VCF only**, **IT only**, **VCF + IT**, **No data**, **Selected**  

These base colors and outlines change automatically by each node’s metadata.  
When hovering over nodes, a **tooltip** appears showing parent/child relations.

---
<h2 id="controls">Pedigree Interaction Controls</h2>

These tools expand, clean, reset, or resize the graph for exploration.

**Expand / Remove / Reset / Wide — core view controls**

- **Expand** — Adds hidden parents/children (+1 generation) to the current node.  
  Useful for exploring connected but unseen relatives.  
  ![expand](/assets/page2_expand.gif)

- **Remove** — Deletes selected node(s) from the current pedigree (base node not removable).  
  ![remove](/assets/page2_remove.gif)

- **Reset** — Returns the pedigree to its initial state.  
  If varieties were added through “Add Option”, the reset also includes those changes.  
  (To fully revert, remove the additions first.)  
  ![reset](/assets/page2_reset.gif)

- **Wide** — Expands the canvas horizontally, enabling easier navigation when many nodes exist.  
  ![wide](/assets/page2_wide.gif)

---
<h2 id="add-option">Add Option Overview</h2>

This section handles adding new resources into the pedigree network,  
combining **pedigree search**, **phenotype search**, and **VCF-based search**.

Note: Phenotype (Page 1) and VCF searches may include varieties **without pedigree**.  
Such nodes appear as **rectangles**, while pedigree-linked nodes remain **circles**.

- **Add** — Opens the **Add Variety** dialog.  
  Search the **pedigree DB** and expand the **±2-generation pedigree** around each selected variety.  
  ![add](/assets/page2_add.gif)

- **Add2 (Phenotype-filtered Add)** — Use Page 1–style **multi-filter for phenotypes** to curate candidates.  
  For each selected variety **with pedigree records**, expand its **±2-generation pedigree**.  
  If no pedigree exists, it remains a **stand-alone candidate**.  
  ![add2](/assets/page2_add2.gif)

- **Add3 (VCF-based Add)** — Pick candidates from the **VCF table**; filter by **ecotype** and **variety group**.  
  For each selected variety **with pedigree records**, expand its **±2-generation pedigree**.  
  If no pedigree exists, it remains a **stand-alone candidate**.  
  ![add3](/assets/page2_add3.gif)

- **Apply** — Commits all queued additions (Add / Add2 / Add3) and **updates the layout**,  
  materializing each selected variety’s **±2-generation pedigree** and **merging overlaps** where lineages intersect.  
  ![apply2](/assets/page2_apply2.gif)

- **Highlight** — Toggles **pedigree-group highlighting** (e.g., *pedigree1*, *pedigree2*).  
  Click a group label to **emphasize only the edges** connecting that pedigree’s members.  
  When expansions cause two pedigrees to intersect, they **automatically merge** (e.g., *pedigree2* merges into *pedigree1*).  
  ![highlight](/assets/page2_highlight.gif)

---

<h2 id="phenotype"> How to Use Phenotype Tab</h2>

Connects the **pedigree (left)** with the **phenotype view (right)**.  
From the pedigree, pick varieties and see **where they fall in the full trait distribution**.  
Only nodes with phenotype data (**IT only** or **VCF + IT**) are applied.

- **Select varieties on the pedigree** — Click one or more nodes that have phenotype data  
  (**IT only / VCF + IT**). The selection is sent to the Phenotype tab.  
  ![pheno_select](/assets/page2_pheno_selet.gif)

- **Open / focus the Phenotype tab (right panel)** — The Phenotype panel lives on the **right**.  
  Because the layout is wide, the UI provides a **switch** to focus that panel.  
  There you can view the **full distribution** and **summary stats** for the chosen trait(s) and selected varieties.  
  ![view_pheno_tab](/assets/page2_view_phenotype_tab.gif)

- **Use the Phenotype tab** — Choose a **trait** to define a color mapping by value range.  
  The app assigns colors to the **selected varieties** accordingly and **propagates the same mapping back to the pedigree** (left),  
  so nodes reflect their trait level in context.  
  ![using_pheno_tab](/assets/page2_using_phenotype_tab.gif)

*Note:* Nodes without phenotype data are skipped

---
<h2 id="gwas">How to Use GWAS Tab</h2>

Available **only for varieties with VCF**.  
Provides **Filter options** (subtrait, MAF, SNP presence, unique) and **Group options** (trait / variant / sample),  
and synchronizes a **bar plot** with **scatter or table** views.

- **Enable GWAS for VCF-backed varieties** — Select one or more nodes that have VCF.  
  The GWAS panel activates and loads trait-linked SNP metadata for those selections.  
  ![gwas_overview](/assets/page2_gwas_overview.gif)

- **Set Filter options** — Adjust **subtrait**, **MAF threshold**, **SNP presence**, and **Unique**.  
  The **bar plot** updates to show the **count of traits** satisfying the current filters.

  **Definitions**
  - **SNP presence**: *presence ≥ N samples*. If multiple samples are selected, a variant passes when it satisfies the **smallest N** among the selected samples.  
  - **Unique**: returns variants **specific to each sample** (sample-unique variants).

- **Apply Group options** — Choose grouping **mode** (trait / variant / sample), define groups, then apply.  
  The **bar plot** reflects **per-group trait counts** that meet the active conditions.

- **Switch View: Bar ↔ Scatter/Table** — Toggle to **GWAS scatter** or **table**.  
  These views inherit the current **filters/groups** and let you inspect SNP positions and summaries.

- **Capture variants & reflect to Pedigree** — Click a **node** in the scatter or a **row** in the table to collect a **Variant ID list**.  
  The selected variants can be **sent back to the pedigree** (overlay/highlight) for lineage-aware interpretation.