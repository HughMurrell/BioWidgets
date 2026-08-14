# aliViz User Guide

aliViz is a browser-based alignment and phylogeny viewer. It supports loading alignments, grouping and founders, inferring or importing trees, clustering, epitope definition and logos, and 3D structure viewing.

Toolbar controls below appear in the **same left-to-right order** as in the app.

---

## Shared concepts

These values appear in the Phylogeny header and in several tools:

- **target-dpi** — Days post infection used for **etd**, k-DPI Auto stopping, and SVG **cluster target bars**. Defaults from the alignment **filename** (`dpi+N` / `dpi-N`, case-insensitive; largest `N` if several; default **14** if none). After a successful **k-DPIs → Accept**, the dialog’s **Target DPI** is remembered for subsequent SVG plots (and shown in the Phylogeny header). Also updates on prune, NT→AA conversion, and when you edit DPI days in the SVG export dialog (until the next k-DPI Accept overrides it).
- **dsr** — Daily substitution rate (default **6.5×10⁻⁵** substitutions/site/day). Resets on page reload. Set with the **dsr** button. Used for tree-based DPI, **etd**, and SVG target-bar length **τ = target-dpi × dsr**.
- **etd** — Expected mean pairwise tip distance under an early-infection / strict-clock star model: **etd = 2 × dsr × target-dpi**. Used as the long-branch noise threshold in k-DPIs (incoming branch **> etd**).
- **Tree-based DPI** — For a set of tips: **dpi = mean pairwise tree path distance / (2 × dsr)** (rounded to nearest day for labels). **Alignment sequence distances are never used for DPI.**
- **Default dsr** — The built-in **6.5×10⁻⁵** comes from early-infection parameters assembled by Keele *et al.* (2008): reverse-transcriptase error rate **μ = 2.16×10⁻⁵** per site per replication, generation time **t_gen = 2** days, and basic reproductive ratio **R₀ = 6**, via

  **dsr = (μ × R₀) / t_gen = (2.16×10⁻⁵ × 6) / 2 ≈ 6.5×10⁻⁵**.

  Here **μ / t_gen** is the baseline lineage clock (one replication per generation), and the factor **R₀** treats Keele’s expanding early-infection regime as raising the *effective* daily substitution rate used when converting tree path length into days (**dpî = δ̄ / (2 × dsr)**, **etd = 2 × dsr × DPI**). Override with the **dsr** button if a different calibration is preferred.

**Functionality tags** (set on every load):

- **AA:** ungapped sequence starts with **M**, ends with **`*`**, no internal **`*`**.
- **NT:** translate frame 1, then apply the AA rules.
- Load summary above the filename: **`s=…, f=…, nf=…`** (all sequences). SVG title lines use **`s,f,nf` for tree leaves only** (REF, subtype, PDB excluded).
- Tree tips are cluster-coloured diamonds; **non-functional** tips get a **red oval**. Red is reserved for that marker (not used in the group/cluster palette).
- On exported **linear** SVGs with clusters: grey dots at ordinary internal nodes; larger filled dots at each cluster’s LCA (cluster colour); solid vertical **target bars** at **LCA depth + τ** (**τ = target-dpi × dsr**), spanning that cluster’s tip rows.
---

## 1. Info

### ?
- Opens a short help overlay. Link to this full guide is in the overlay footer.

---

## 2. Sequences

### Choose file
- Load a **FASTA** or **FASTQ** alignment. The filename appears beside the button.
- **Character sanitization:** Keep **A–Z**, **`-`**, **`*`**; any other character → **`X`** (with an alert of replacements).

### Has reference
- **Default: checked.** Sequence 1 is **Reference (REF)**.
- **Unchecked:** Choose file opens a **dual-load** dialog (alignment + optional separate reference FASTA). The reference is prepended for display and analysis.

### Has subtype
- **Default: unchecked.** When checked, the sequence after REF (or sequence 1 if no REF) is **SubType**, labelled **[Subtype]**, and excluded from tree inference unless it is also the founder.
- Toggling after load **clears grouping, tree, and clustering**.

### Collapse
- Read an integer **collapse count** from a field in each sequence name (delimiter + 1-based field number). Counts weight consensus, medoid, and related tallies.
- Leave the field number empty to clear collapse.

### Group
- Split each name by a **delimiter** and take a **1-based field** as the group label. REF, SubType (if any), and PDB chains are excluded.
- Legend: **`label (n)`** until a phylogeny exists; then **`label (n, dpi)`** with tree-based DPI when available.

### Sort
- Reorder by current grouping (or by name if ungrouped). REF and SubType stay first.

### Prune
- Remove sequences belonging to selected groups. Remaining group IDs are renumbered. The tree is cleared; group legend returns to counts only.

### Add Founder
- Requires grouping (for group-based founders).
- **Founder source (default: Medoid):**
  - **Consensus of a group** — New majority-per-column sequence named `consensus_of_<label>` (or `consensus_of_subtype`).
  - **Medoid sequence of a group** — Existing sequence minimizing total padded Hamming distance to others in the group; keeps its name; shown as **[Founder]**.
- Clears any existing tree. **Infer** enables once a founder is set. The founder is not forced into cluster 0; cluster 0 is reserved for subtype when present.

---

## 3. Phylogeny

Header shows **`dsr=…, target-dpi=…`**.

### dsr
- Set the session daily substitution rate. **Use default** restores **6.5×10⁻⁵** (Keele *et al.* 2008; see Shared concepts). Applying refreshes tree-based group/cluster DPI labels when a tree is present, and reseeds SVG Mutations/day from the new dsr.
### Infer
- **Requires a founder** (button stays disabled until then).
- Builds or loads a phylogeny for sample sequences (excludes REF, subtype when present, and PDB chains; includes founder).
- **Methods:**
  - **Load inferred tree (Newick)** (dialog default) — `.nwk`, `.newick`, `.tree`, `.tre`, `.treefile`, `.txt`. Leaf names must match the alignment after normalization (`_cl-*` ignored).
  - **FastTree (bioWASM)** — NT/AA model and speed options.
  - **Neighbor-Joining (Internal)** — PhyloTools pairwise distances → NJ.
- **Root & ladderize** (checked by default): after infer/load, reroot on the **founder** then ladderize **by depth**. Uncheck to do those steps manually.
- After a tree exists: if sequences are grouped, legend gains tree-based **dpi** per group.

### Reroot
- Reroot on **[SUBTYPE]** or **[Founder]** (dialog). Topology and branch lengths are preserved; only the root moves.

### Ladderize
- Reorder tree children and the sequence list to match leaf order.
- **By Depth** (default): sort children by max root-to-leaf depth in the subtree.
- **By Weight:** sort by leaf count.
- Enables Histogram, Cluster, and tree downloads after ladderizing.

### Histogram
- Histogram of **root-to-leaf** branch-length distances. Requires a ladderized tree.

### Cluster
Opens a method selector, then the chosen interface. Sample sequences only (REF / subtype / PDB excluded).

**Methods:** **None**, **Hierarchical (tree-clade)**, **Hierarchical (tree-cut)**, **k-DPIs (divisive max-depth)**, **UMAP**, **MDS** (default).

#### None
- Clears clustering, strips `_cl-*` from names and tree nodes. Group colours remain.

#### Auto (hierarchical / MDS / UMAP)
- Dialogs with an **Auto** button search discrete slider ticks and maximize the selected index:
  - **Calinski–Harabasz** (default) or **Ball-Hall-Adapted** (radio).
- Tree-clade: every branch-length threshold tick.
- Tree-cut: all depth triples with d1 ≤ d2 ≤ d3.
- MDS/UMAP: every **eps** tick for the current min neighbors.

#### Hierarchical (tree-clade)
- New clade when incoming **branch length** exceeds a threshold (DFS). **Min leaves** → noise (`_cl-n`). Accept writes `_cl-*`, tree-based cluster DPI, legend and tip colours. Cancel strips clustering and keeps groups.

#### Hierarchical (tree-cut)
- Three **depth** cutoffs (Level 1 ≤ 2 ≤ 3). Leaves assigned by crossing depth lines (BFS). Min leaves → noise. Accept / Cancel as above.

#### k-DPIs (divisive max-depth)
- Repeatedly split clusters on the **phylogeny**. Each candidate cut uses an **internal node**: **right** = descendant tips in the cluster, **left** = remainder. Choose the cut that minimises **max(depth₁, depth₂)** (LCA-subtree height). Only cuts where **both** children have size ≥ **min cluster size** are considered (a peel of a single tip therefore requires leaf size ≥ **2 × min size**).
- **Target DPI** (dialog field; default from filename `dpi±N` or **14**): **Auto** only splits leaves with **dpi > Target DPI** and stops when every remaining cluster has **dpi ≤ Target DPI** (or **Max splits** is reached). Changing Target DPI does **not** recluster by itself — use **Reset**, then **Auto** (or **+**) with the new target, then **Accept**.
- **+** / **−** change k; **Max splits** (default **9**) caps **+** and **Auto**. Both apply the best depth bipartition when one exists. **+** picks the non-noise leaf with largest **DPI** that can still yield two min-sized children.
- **Min cluster size** (default **1**): after splitting, leaves with **n < min size** become **noise** (`_cl-n`) and are never split further. Singletons are therefore noise when min size &gt; 1; with the default of 1 they remain ordinary clusters.
- **Remove long-branch tips as noise** (default on): at open / Reset / Auto start, tips with incoming branch **> etd** (**etd = 2 × dsr × Target DPI**) become pre-noise and are excluded from splitting.
- Bubble tree: **root left**, remainder **up**, descendants **down** (like the linear phylogeny); leaf label **`n, dpi`**; radius ∝ log(n).
- **Accept:** writes `_cl-<id>` / `_cl-n`, legend **`id (n, dpi)`** and **Noise (n)**, and **remembers Target DPI** for subsequent SVG cluster target bars (Phylogeny header **target-dpi** updates). **Cancel** drops the preview without changing the remembered target.
#### UMAP
- UMAP on tree pairwise distances → 2D, then **DBSCAN**. Parameters: nNeighbors, Spread, Min Distance, plus eps / Min Neighbors. **Estimate DPI** and **Apply** as for MDS.

#### MDS
- Classical MDS of tree pairwise distances → 2D (`B = −0.5 · H · D² · H`), then **DBSCAN**. **Estimate DPI** shows tree-based **`clusterId (n, dpi)`**. **Apply** stamps clusters and legend DPI; **Close** discards.

#### DBSCAN (MDS / UMAP)
- Neighbors within **eps**; cores need ≥ **minPts** neighbors; else noise (−1).

#### Calinski–Harabasz (tree-based)
- Scores clusterings for Auto. Between / within terms use MRCA and root distances on the tree; result divided by **2<sup>k</sup>**.

#### Ball–Hall-Adapted (tree-based)
- **BH<sub>adapted</sub> = [ BH(1) / BH(k) ] / 2<sup>k</sup>**, with BH = sum of squared leaf–MRCA distances. Shown beside CH in hierarchical dialogs.

---

## 4. Download

### fasta
- Download the current view (NT or AA). Filename: **`<alignment>_gr-<n>_cl-<n>.fasta`** (clusters exclude noise).

### tree
- Export Newick. Filename: **`<alignment>_gr-<n>_cl-<n>.nwk`**. Requires a ladderized tree.

### svg
- Export a vector tree figure. Same filename stem with **`.svg`**. Dialog: **Geometry** (linear / circular) and **Scale**, then Export.

**Geometry**
- **Linear** — Rectangular cladogram.
- **Circular** — Radial orthogonal layout (half scale per unit vs linear).

**Scale** (branch region ≈ 520 px; 1/16 mark at 32.5 px)
- **Auto** — Fit the tree to the plot. If the filename contains `dpi±N`, also draws a DPI scale annotation and fits so that mark (and any cluster target bars) stay in view.
- **Fixed** — Branch length per 1 cm (dropdown). Overflow aborts with a suggestion to pick a larger value or Auto.
- **DPI** — Place expected depth **DPI × Mutations/day** at the 1/16 mark. DPI days default from filename (or last k-DPI Accept when that has updated the SVG DPI field); Mutations/day seeded from session **dsr**. Never aborts on overflow (clips at the right border; linear clipped tips use a **right-half diamond**; non-functional tips keep the red oval). Gray dotted line at expected depth.

**SVG extras**
- Title: filename; second line = infer method, leaf **`s,f,nf`**, max depth, last cluster method.
- Scale bar on the tree panel (fixed: 1 cm; DPI: 32.5 px linear / 16.25 px circular).
- Legend panel: groups and clusters (numeric cluster IDs).
- **Linear + clusters:** per-cluster solid target bars at **LCA + τ** with **τ = target-dpi × dsr** (target-dpi from last **k-DPI Accept** when set, otherwise filename / SVG DPI field). Bars use the cluster colour and span only that cluster’s tip rows. Small grey markers on ordinary internal nodes; larger cluster-coloured markers on each cluster LCA.
---

## 5. Reset

### Clear Alignment
- Remove all sequences and reset related state (groups, clusters, tree, founder, etc.).

### Clear Tree
- Remove the phylogeny and clustering. Disables Histogram, Cluster, and tree downloads.

---

## 6. View Mode

### NT / AA
- Switch nucleotide vs amino-acid display. AA translates with the standard genetic code (frame 1).

### AA palette
- **Alignment (IUPAC)** (default) or **Biochemical**.

### Highlighter
- Highlight mismatches vs **Off**, **Founder**, **SubType**, or **Reference**.

---

## 7. Epitopes

### Load Epitopes
- CSV: first column = name; further columns = regions `start:end` or a single position. Example:

  `VRC01,197:198,230,276,278:282,365:371,427:428,430,455:463,465,469,471:474`

- Merges with existing epitopes (same name overwrites). Loading in NT mode switches the view to AA.

### Select Epitope
- Restrict the alignment display to one epitope’s regions, or **None** for all columns.

### Show Logo
- Sequence logo for the selected epitope (requires a selection).

### Export Epitopes
- Download current epitopes as CSV (same format as load). Disabled when empty.

### New epitope (from View 3D)
- Distance-based epitope: residues on an **intrinsic** chain within a cutoff (Å) of an **extrinsic** chain. Name `<pdbCode>_<cutoff>A`. Mapped to alignment columns via the reference.

---

## 8. 3D Structure

### Load file
- Load PDB/CIF from disk, or enter a 4-character PDB ID and **Fetch** from RCSB. Chains become sequences (`xxxx_Chain_y [PDB_y]`) aligned to the reference in AA mode.

### View 3D
- 3Dmol.js viewer. Chain categories:
  - **Intrinsic (≥50%)** — In the alignment and ≥50% coverage/identity to REF.
  - **Extrinsic (<50%)** — In the alignment but below that threshold.
  - **Superficial** — Not in the alignment panel.
- Epitope colouring on intrinsic chains; residue click can open a partial logo. **New epitope**, Full screen, Export PNG, Close.

---

## 9. Color legend

- **Groups:** `label (n)` after Group; `label (n, dpi)` after a phylogeny when DPI is available.
- **Clusters:** After Accept/Apply, `clusterId (n, dpi)` when available; **Noise** when present.
- Phylogenetic scale bars appear only on **exported SVG** figures, not in the floating legend.

---

## 10. Algorithms and formulas (quick reference)

| Feature | Summary |
|---------|---------|
| Load sanitize | A–Z, `-`, `*` kept; else → `X`. |
| Functionality | AA: M…`*`, no internal `*`; NT: frame-1 translate then same; `s,f,nf` above filename. |
| DPI / etd | dpi = mean tree path / (2×dsr); etd = 2×dsr×target-dpi; default dsr = 6.5×10⁻⁵ (Keele 2008: μ R₀ / t_gen); dsr via **dsr** button. |
| Group | Name field → group; legend gains tree DPI after phylogeny. |
| Collapse | Integer count from name field; weights consensus/medoid. |
| Founder | Consensus (new row) or medoid (existing name + `[Founder]`). |
| Infer | NJ / FastTree / Load Newick; optional root on founder + ladderize by depth. |
| Reroot | On subtype or founder. |
| Ladderize | By depth or by weight. |
| SVG scale | Auto fit (+ optional filename DPI overlay); Fixed length/cm; DPI at 1/16 mark (clip, half diamond). |
| SVG cluster bars | Linear: LCA + τ, τ = target-dpi × dsr; target-dpi from last k-DPI Accept when set. |
| Tree-clade | Edge length threshold; min leaves → noise. |
| Tree-cut | Three depth cuts; min leaves → noise. |
| k-DPIs | Min max-child depth; Target DPI field; Auto only splits DPI &gt; target; Reset+recluster after changing target; Accept remembers target for SVG bars; min size default 1 → noise; long-branch pre-noise via etd; bubble tree root-left. |
| MDS / UMAP | Tree distances → 2D → DBSCAN; CH / Ball-Hall Auto on eps. |
| CH | Tree MRCA/root form; ÷ 2<sup>k</sup>. |
| Ball-Hall-Adapted | BH(1)/BH(k) ÷ 2<sup>k</sup>. |
| Epitope CSV | name, regions… |
| Distance epitope | Intrinsic–extrinsic min distance &lt; cutoff (Å). |

---

*End of aliViz User Guide.*
