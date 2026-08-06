# aliViz User Guide

aliViz is a bioinformatics alignment and phylogeny viewer. It supports loading alignments, inferring or importing trees, grouping sequences, clustering (tree-based and projection-based), epitope definition and logos, and 3D structure viewing with epitope coloring.

---

## 1. Loading and view controls

### Choose file
- **Function:** Load an alignment from a **FASTA** or **FASTQ** file.
- **Usage:** Click **Choose file**, select your alignment. The loaded filename appears in the text box beside the button (replacing the initial “No file chosen” placeholder).
- **Reference:** If **Has reference** is checked (see below), the **first** sequence is treated as **Reference (REF)**. If unchecked, you load the alignment and reference separately (dual-load dialog).
- **Subtype:** If **Has subtype** is checked (see below), the **second** sequence is treated as **SubType** when a reference is present (otherwise sequence 1). If unchecked, no sequence is designated as subtype unless you add one later by other means.

### Has reference
- **Function:** Tell aliViz whether the alignment file already includes the reference as its **first** sequence.
- **Default:** **Unchecked** (many workflows load the reference from a separate FASTA).
- **When checked:** Choose file loads the alignment directly; sequence 1 is **Reference (REF)**.
- **When unchecked:** Choose file opens a **dual-load** dialog: pick the alignment, then optionally pick a reference FASTA (or reuse a reference loaded earlier in the session). The reference is prepended to the alignment for display and analysis.

### Has subtype
- **Function:** Tell aliViz whether the alignment includes a dedicated subtype sequence in position 2 (after reference).
- **Default:** **Unchecked** (many alignments have no subtype).
- **When checked:** Sequence 2 is the subtype; it receives the **[Subtype]** label and is excluded from tree inference (unless it is also the founder).
- **When unchecked:** No subtype index is assigned; subtype labeling and subtype-specific clustering behaviour are disabled.
- **After load:** You can toggle this checkbox. Changing it **clears grouping, tree inference, and clustering** and resets related state so you can re-run those steps with the new subtype setting.

### Character sanitization (FASTA/FASTQ load)
- **Allowed characters:** Uppercase letters **A–Z**, gap **`-`**, and stop codon **`*`** (common in amino acid alignments).
- **Non-standard symbols:** Any other character in a sequence line is **replaced with `X`** (alignment length is unchanged). You receive an alert listing what was replaced (e.g. `? → X (12)`).
- **Nucleotide alignments:** Standard IUPAC nucleotide letters are within A–Z and are kept as-is.

### Sequence functionality (functional / non-functional)
On every alignment load, aliViz tests each sequence for **functionality** and stores a **`functional`** tag (`true` or `false`) on that sequence.

- **Amino acid (AA) alignment:** A sequence is **functional** when its **ungapped** sequence (gaps ignored):
  - **Starts with M** (start codon),
  - **Ends with `*`** (stop codon),
  - Has **no internal `*`** (no premature stop codons).
- **Nucleotide (NT) alignment:** Each row is translated using **reading frame 1** (same default as at load), then tested with the AA rules above.
- **Non-functional:** Any sequence that fails one or more of these rules.

**Load summary (`s`, `f`, `nf`, `dpi`, `dsr`, `etd`):** Above the alignment filename you will see a compact report, e.g. **`s=56, f=25, nf=31, dpi=24, dsr=0.000079, etd=0.003792`**, meaning 56 sequences total, 25 functional, 31 non-functional, **24** days post infection (DPI), daily substitution rate **dsr**, and **expected mean pairwise tip distance** **etd**.

- **DPI** is the same value used as the default in the tree SVG **DPI** scale dialog (see §3). It is set from the alignment **filename** when the file loads: aliViz looks for tokens like `dpi+24` or `dpi-24` (case-insensitive); if several tokens are present, the **largest** day count is used; if none are found, DPI defaults to **14**. It updates when you **Prune** sequences or when NT rows are converted to AA (e.g. after loading epitopes), and when you change the DPI (days) field in the SVG export dialog.
- **dsr** is the current **daily substitution rate** (default **7.9×10⁻⁵**). It resets to this default whenever the page is reloaded. If the alignment **filename** contains a DPI token, a dialog offers to keep the default or set another **dsr** for the session; you can also change **dsr** later in the Group dialog. Changing **dsr** updates this summary and **etd**.
- **Default dsr estimate:** The built-in default **dsr = 7.9×10⁻⁵** substitutions/site/day was obtained from an **early mutation rate** analysis using the companion script **`EarlyMutationRateEstimate.py`**. That script parses exported aliViz tree SVGs: for each file it reads **DPI** (from the filename), **MaxDepth** (from the SVG title line), and **CL** (cluster count from the filename). It selects records with **CL = 1** and **DPI ≤ 175**, plots **MaxDepth** versus **DPI**, and fits a regression line through the origin. The **slope** of the no-outlier fit (after excluding points with large residuals from an initial through-origin fit) is taken as the estimate of **dsr** (MaxDepth ≈ dsr × DPI under a strict clock). Re-run the script on your own SVG directory to reproduce or update the estimate.
- **etd (expected mean pairwise tip distance):** under an early-infection / strict-clock **star** approximation, mean tip-to-tip path length equals the tree diameter and is about twice the expected root→tip depth, so  
  **etd = 2 × dsr × DPI**  
  (e.g. DPI = 24 and dsr = 7.9×10⁻⁵ → etd = 0.003792 substitutions/site). This is an expected scale for the phylogeny, not a measured tree statistic.

The summary resets to **`s=0, f=0, nf=0, dpi=14, dsr=0.000079, etd=…`** when the alignment is cleared (or on page reload before a file is loaded).

**Tree indicators (interactive panel and SVG export):** Whenever a tree is drawn— in the **tree panel** beside the alignment, or in an exported **SVG** figure— leaf tips are shown as **diamonds** coloured by **cluster** (when clustering is active). **Non-functional** leaves are marked with a **red oval** (`#dc2626`) drawn around the diamond. The oval appears on:
- Full tip diamonds (normal tips),
- **Right-half diamonds** on linear SVG exports when a tip is **clipped** at the plot boundary (DPI scale overflow).

Functional leaves show only the cluster-coloured diamond with no oval. **Red is reserved** for this non-functional marker: the shared group/cluster colour palette does **not** use red or pink (palette indices 5 and 10 are sky and slate instead), so the oval remains visually distinct from cluster colours.

**SVG title line:** Exported tree SVG figures include **`s=…, f=…, nf=…`** counts on the **second title line**, placed after the **Infer** method and before **Max depth**. These counts reflect **visible tree leaves only** (reference, subtype, and PDB chains are excluded because they are not plotted on the tree), e.g. `Infer: FastTree (bioWASM) | s=54, f=24, nf=30 | Max depth: 0.014 | Cluster: …`. The load summary above the filename still counts **all** sequences in the alignment file.

### View Mode: NT / AA
- **Function:** Switch between **Nucleotide (NT)** and **Amino acid (AA)** view.
- **Algorithm:** In AA mode, nucleotides are translated using the standard genetic code (frame can be set if frame selector is enabled). Gaps and invalid codons produce gap characters or `X`.
- **Modes:** **NT** (default), **AA**.

### AA palette
- **Function:** Choose the amino acid color scheme.
- **Options:** *Alignment (IUPAC)* (default) or *Biochemical* (e.g. by property).

### Highlighter
- **Function:** Highlight differences in the alignment with respect to a chosen reference.
- **Options:** **Off** (default), **Founder**, **SubType**, **Reference**. Mismatches against the selected sequence are highlighted.
- **Founder:** Works for both a **consensus** founder (`consensus_of_…` sequence) and a **medoid** founder (an existing sequence marked as founder; see §2).

---

## 2. Grouping and sequences

### Group
- **Function:** Assign sequences to groups using a delimiter and a field in the sequence name.
- **Algorithm:** You specify a **delimiter** (e.g. `_`) and a **field number** (1-based). Each sequence name is split by the delimiter; the value at that field becomes the group label. REF, SubType (if present), and PDB chain sequences are excluded from grouping. Groups are assigned unique IDs and used for coloring and legend.
- **DPI estimate (per group):** After grouping, aliViz estimates **days post infection** for each group from the **mean** pairwise **p-distance** (**api**: average fraction of differing non-gap sites among all pairs in the group) and the **daily substitution rate** **dsr** from the Group dialog (default **7.9×10⁻⁵** for the session; not persisted across page reloads). Under a Poisson substitution process, the linear rule `dpi = api / (2 × dsr)` is correct only in the **infinite-sites** limit (no multiple hits). When sites can change more than once, observed `api` saturates and underestimates divergence. aliViz therefore applies a **Jukes–Cantor (JC69)** multiple-hit correction, then converts corrected divergence to time with the same clock:
  - **d** = −(3/4) ln(1 − 4·api/3) (undefined when api ≥ 0.75)
  - **dpi** = d / (2 × dsr), rounded to the nearest day
  - Shown in the legend as **`label (n, dpi)`**, e.g. **`2000 (12, 56)`** (n = sequence count). Groups with fewer than two sequences have no DPI estimate. Estimates are recomputed after **Prune**.
- **Result:** `state.sequenceGroups` (name → group ID), `state.groupNames` (group ID → label), and `state.groupDpiEstimates` (group ID → dpi days) are set; sequence names are colored by group in the name panel.

### Sort
- **Function:** Reorder sequences by current grouping (or by name if no grouping), keeping REF first and SubType second (when a subtype exists).
- **Algorithm:** REF and SubType stay at the top; other sequences are sorted by group ID (from `state.sequenceGroups`), then by name within each group. PDB chains and founder are handled in the sort order.

### Prune
- **Function:** Remove sequences that belong to selected groups.
- **Usage:** Open the Prune overlay, select groups to **remove**, then Apply. Removed sequences are deleted from the alignment. Group IDs are renumbered to 0, 1, 2, … for the remaining groups.

### Add Founder
- **Function:** Designate a **founder** sequence for a chosen group (or from the subtype sequence when offered in the group list).
- **Founder source (default: Medoid):**
  - **Consensus of a group:** Builds a **new** consensus sequence (majority character per column, including gaps) from all sequences in the group (excluding REF and subtype). It is inserted into the alignment with a name such as `consensus_of_<groupLabel>` (or `consensus_of_subtype` when using the subtype option in consensus mode).
  - **Medoid sequence of a group (default):** Selects an **existing** sequence in the group that minimizes total **padded Hamming distance** to all other sequences in that group. The sequence **keeps its original name**; it is marked as founder and shown with the **[Founder]** label in the name column and tree.
- **Switching modes:** Replacing a consensus founder removes or overwrites the synthetic `consensus_of_` row; switching to medoid on a group removes a previous synthetic founder if present and tags the medoid sequence instead.
- **Clustering:** The founder is **not** forced into cluster 0. When **no subtype** is present, cluster IDs start at **1** (cluster 0 is reserved for subtype when it exists). A medoid founder participates in clustering like any other sample sequence; the **[Founder]** label is preserved after clustering (name matching tolerates `_cl-*` suffixes added by clustering).
- **Usage:** Group sequences first (for group-based founders). Open **Add Founder**, choose **Founder source**, select a group (or **[SubType]** when available), then **Add Founder**. Any existing tree is cleared; **Infer** becomes enabled once the founder is settled.

---

## 3. Tree inference and manipulation

### Infer
- **Requirement:** A **founder** must be designated first (**Add Founder**). The **Infer** button stays disabled until then.
- **Function:** Build or load a phylogeny for the current alignment (excluding REF, subtype when present, and PDB chains; **including** founder sequences, including medoid founders).
- **Methods:**
  - **Load inferred tree (Newick file)** (default in the infer dialog): Opens a file picker for an existing **Newick** tree. Accepted extensions: **`.nwk`**, **`.newick`**, **`.tree`**, **`.treefile`**, **`.txt`**.
  - **FastTree**: Uses FastTree (via bioWASM) on the alignment. Model options appear for NT and AA.
  - **Neighbor Joining (NJ):** Uses PhyloTools (pairwise distances → NJ tree).
- **Root & ladderize:** The infer dialog includes a **Root & ladderize** checkbox (on by default). When checked, after the tree is inferred or loaded aliViz automatically **reroots on the founder** and then **ladderizes by depth** (same as the Reroot and Ladderize tools). Uncheck to leave the raw inferred/loaded tree and run those steps manually.
- **Loading an external tree:** Leaf names in the Newick file must match alignment sequence names after normalization (cluster suffixes `_cl-*` are ignored for comparison). **Reference**, **subtype** (if present), and **PDB chains** are excluded from the required leaf set. A medoid founder keeps its sample name in the alignment, so tree leaves should use that name—not a separate `consensus_of_` name. If names do not match, loading is aborted with a message listing missing or extra leaves.
- **Mean pairwise distance:** After a tree is inferred or loaded, aliViz computes the **mean pairwise tip-to-tip path length** (branch-length units) and shows it beside the **Phylogeny** control label as **`mean=…`**. It clears when the tree is removed.

### Reroot
- **Function:** Reroot the tree on the **Founder** sequence.
- **Targets:** **Founder** (requires a founder to be defined—consensus or medoid).
- **Algorithm:** The founder leaf is located and the tree is rerooted at the edge leading to it so that the founder is the outgroup. Branch lengths and topology are preserved; only the root position changes.

### Ladderize
- **Function:** Reorder the tree and the sequence list so that the tree is ladderized and sequences match leaf order.
- **Modes:**
  - **By Depth** (default): At each node, compute the maximum root-to-leaf distance in each child’s subtree (using branch lengths). Sort children by that depth (ascending). Shallower subtrees appear first.
  - **By Weight:** At each node, sort children by number of leaves (ascending). Lighter subtrees appear first.
- **Result:** Tree drawing order and `state.viewSequences` are updated to follow the new leaf order. Histogram, Cluster, and tree export (Newick **tree** and **svg**) are enabled after ladderizing.

### Histogram
- **Function:** Show a histogram of **root-to-leaf distances** (sum of branch lengths from root to each leaf).
- **Algorithm:** For each leaf, the path from root to leaf is traversed and branch lengths (`node.len`) are summed. Distances are binned (number of bins between 10 and 30, about √n). Bars are drawn for each bin count.
- **Requirement:** Tree must be loaded and ladderized.

### tree (Export tree)
- **Function:** Export the current tree in **Newick** format.
- **Usage:** Click “tree” to download a Newick file containing the current tree with names and branch lengths (if present).

### svg (Export tree as SVG)
- **Function:** Download a **vector (SVG)** figure of the current tree for publications or slides.
- **Requirement:** A tree must be loaded and **ladderized** (same as the Newick export).
- **Usage:** Click **svg** in the phylogeny toolbar. A dialog offers **Geometry** and **Scale**, then **Export** to save `tree.svg`. **Cancel** closes without downloading.

#### Geometry
- **Linear:** Rectangular cladogram (horizontal branches, vertical backbone at internal nodes, diamond tips, dashed connectors to labels).
- **Circular:** Radial layout with **orthogonal** branches (circumferential arc + radial segment), outward-pointing tip diamonds, and dotted connectors to radial labels.

#### Scale
The plot has a **520 px** branch region (≈ **13.76 cm** at 96 px/inch) available for branch lengths; its 1/16 mark is at **32.5 px**.

- **Auto (default):** Scales the tree to fit the plot area (linear branch span ≈ 520 px; circular outer branch radius ≈ 260 px, with branch lengths in circular mode drawn at **half** the linear scale per unit).
- **Fixed:** Uses a chosen **branch length per 1 cm** of plot (dropdown: 0.5, 0.1, 0.05, 0.01, 0.005, 0.001, 0.0005, 0.0001; default **0.001**) so trees from different datasets are comparable. The tree is **not** shrunk to fit; shallow trees use only part of the space. If the tree would extend beyond the fixed plot width, export **aborts** with a message suggesting a **larger** value in the dropdown (more compact drawing) or **Auto** scale.
- **DPI (days post infection):** Scales so an **expected** branch length lands at the **1/16 mark (32.5 px)** of the branch region. You supply two parameters:
  - **DPI (days):** number of days post infection. The default is taken from the loaded alignment **filename**: tokens matching `dpi+N` or `dpi-N` (case-insensitive) are parsed, and if several are present the **largest** `N` is used; if the filename has no such token, the default is **14**. The load summary’s `dpi=` value and this dialog field stay in sync. You can still edit the field before export.
  - **Mutations/day (MPD):** expected substitution rate per day. Seeded from the session **dsr** at alignment load (including after the on-load dsr dialog); default **7.9×10⁻⁵**. You can still edit it in this dialog before export.
  - The expected branch length is **DPI × MPD** (e.g. 14 × 7.9×10⁻⁵ = 0.001106), and the scale is set so this value maps to 32.5 px (i.e. `px per unit = 32.5 / (DPI × MPD)`). In circular mode the same value maps to a **16.25 px** radius (half, matching the Fixed circular convention).
  - **No overflow abort:** unlike Fixed scale, DPI **never** aborts. Branches longer than the plot width are truncated at the **right border** (so they cannot run into the sequence-name column); nothing is clipped on the left. In circular mode over-long branches are clamped to the outer radius.
  - **Clipped-branch marker (half diamond):** on **linear** SVG exports, any horizontal branch that would extend past the right border is drawn only up to that boundary. If the tip falls past the boundary, the full diamond is replaced by the **right half** of the diamond (apex on the boundary, pointing right), filled with the same **cluster** color as an unclipped tip. **Non-functional** clipped tips also receive the **red oval** around the half diamond. The dashed connector to the sequence name still starts from the boundary.
  - **Expected-depth marker:** a **gray vertical dotted line** is drawn across the plot at the expected depth (the end of the scale bar, 32.5 px). Tips to the **right** are mutating **faster** than expected; tips to the **left** are mutating **slower** than expected.

#### Titles, scale bar, and legend (SVG only)
- **Title (two lines):** Full alignment **filename** (no truncation); second line includes **inference method**, **`s=…, f=…, nf=…`** (functionality counts for **tree leaves only**—reference, subtype, and PDB chains excluded), **max tree depth** (3 significant figures), and **last clustering method** (or “No clustering”). Linear layout: titles are **left-aligned** with the scale bar; circular: titles are **centred**.
- **Phylogenetic scale bar:** Drawn on the **tree** panel (top-left for linear, top-right for circular—not in the floating HTML legend). Shows branch-length units; in **fixed** mode the bar is **1 cm** long with a numeric label matching the selected scale (circular bar length is **half** the linear bar for the same scale value). In **DPI** mode the bar is **32.5 px** long (linear; 16.25 px circular) and is annotated with the parameters and expected depth, e.g. **`(14 DPI @ 0.000079 = 0.001106)`**, updating to match whatever DPI and MPD you enter.
- **Legend panel:** Header **Legend** (not “Color Legend”). **Groups** and **Clusters** as in the app; cluster rows show the **number only** (e.g. `1`, not `Cluster 1`). Page width split **7/8** tree, **1/8** legend.

#### Layout details (unchanged behaviour, for reference)
- Tip markers use **cluster** colors where clustering is active; **non-functional** tips add a **red oval** around the diamond (see §1). Label and connector colors use **group** / special-sequence rules (magenta for REF, subtype, founder, PDB as in the app).
- Circular plot height grows as needed for long radial labels.

### Clear Tree
- **Function:** Remove the current tree. Clustering state is cleared. Histogram, Cluster, and related buttons are disabled.

---

## 4. Clustering (Cluster button)

Clicking **Cluster** opens a method selector, then the corresponding clustering interface.
Each method clusters sample sequences; **Reference**, **Subtype** (if present), and **PDB chains** are excluded from clustering assignments.
- **Methods (selector):** **None**, **Hierarchical (tree-clade clustering)**, **Hierarchical (tree-cut clustering)**, **k-DPIs (divisive max-depth clustering)**, **UMAP**, **MDS** (default).

### None
- **Function:** Remove active clustering.
- **Effect:** Clears cluster state, strips `_cl-*` suffixes from sequence and tree node names, and updates the legend. **Group** colors and group membership are preserved.

### Auto button behaviour (optimizes selected index)
- **Radio buttons (Calinski vs Ball-Hall-Adapted):** choose which index the dialog’s **Auto** button maximizes while it searches slider values. This does not change the clustering method itself; it only changes the scoring metric used by Auto (the radio choice is read when you press **Auto**). Default is **Calinski–Harabasz**; switch to **Ball-Hall-Adapted** to optimize that index instead.
- **Tree-clade Auto:** searches every discrete branch-length threshold slider tick (step = max/50) to maximize the **selected index** for the current **min leaves per cluster**, then sets the slider to the best threshold.
- **Tree-cut Auto:** full grid search over depth triples (d1, d2, d3) with d1 ≤ d2 ≤ d3 on the discrete depth slider ticks (step = max/50), and sets the sliders to the triple that maximizes the **selected index**.
- **MDS/UMAP eps Auto:** evaluates every discrete **radius (eps)** slider tick from min to max (same snapping as the HTML range control; step = max/50) and selects the eps that maximizes the **selected index** (Calinski–Harabasz or Ball-Hall-Adapted radio) for the current **min neighbors**.

### 4.1 Hierarchical (tree-clade clustering)
- **Function:** Define clusters by **clades** on the tree: a new clade (and thus a new cluster) starts when the **incoming branch length** to a node exceeds a threshold.
- **Parameters:**
  - **Branch length threshold:** Minimum branch length that starts a new clade (slider from step to max; step = max/50). Leaves are assigned to the clade that contains them.
  - **Min leaves per cluster:** Clusters with fewer than this many leaves are relabelled as noise (`_cl-na`).
- **Algorithm:** DFS from root. When traversing an edge longer than the threshold, increment cluster ID. Assign each leaf to the current cluster. Subtype (if present) is cluster 0; without subtype, numbering starts at **1**. Then apply min-leaves filter, mark small clusters as noise, and renumber clusters 1, 2, 3, …
- **Auto:** Searches every discrete threshold tick (step to max) to **maximize the selected index** (see §4.6) for the current min-leaves. Sets the slider to the best value.
- **Accept:** Applies the clustering: adds `_cl-<id>` or `_cl-na` to sequence and tree node names, computes cluster DPI estimates, updates the legend and tree colors.
- **Cancel:** Removes all clustering: clears `state.leafClusters`, strips `_cl-*` from names and tree nodes, updates legend and redraws. **Group colors** on sequence names are preserved (group mapping is rekeyed to the stripped names).

### 4.2 Hierarchical (tree-cut clustering)
- **Function:** Cut the tree with **three depth lines** (root-to-node distance). Each region between/above lines defines clusters.
- **Parameters:**
  - **Level 1, 2, 3 Depth:** Sliders (0 to max root-to-leaf distance; step = max/50). Order is enforced: Level 1 ≤ Level 2 ≤ Level 3.
  - **Min leaves per cluster:** As in tree-clade; clusters below this size become noise.
- **Algorithm:** Leaves with depth &lt; Level 1 = cluster 1. BFS from root: when a node’s depth crosses Level 1 (or 2 or 3), all leaves in its subtree with depth ≥ that level are assigned the next cluster ID. Then min-leaves → noise, renumber clusters.
- **Auto:** **Full grid search** over all triples (d1, d2, d3) with d1 ≤ d2 ≤ d3 on the discrete depth slider ticks. For each triple, computes the cluster map and the selected index; chooses the triple that maximizes it and sets the three sliders.
- **Accept / Cancel:** Same idea as tree-clade: Accept writes cluster tags, computes cluster DPI estimates, and updates legend/tree; Cancel clears clustering and strips `_cl-*` while keeping group colors.

### 4.3 k-DPIs (divisive max-depth clustering)
- **Function:** Partition sample sequences by repeatedly splitting clusters on the **phylogeny**. Each split cuts at an **internal tree node** and keeps the bipartition that minimises **max(depth₁, depth₂)**, where **depth** is the height of a tip set from its **LCA** (max path length LCA→tip). Bubble / legend **DPI** labels are separate from that depth score: they use **mean pairwise tree path distance**, **dpi = d / (2 × dsr)** (no JC69; same distances as the CSV downloaded on open). Group / MDS “Estimate DPI” still use sequence p-distance + JC69.
- **Requires:** A loaded or inferred **tree** whose leaves match the sample sequences.
- **Parameters:**
  - **k** starts at **1**. Use **+** / **−** to increase or decrease k. Each **+** splits the current leaf with the **largest max depth** (among leaves large enough to split); each **−** undoes the last split. Split history is remembered so decreases restore the prior partition. **+** and **Auto** stop after **Max splits** (default **9**).
  - **Max splits** (default **9**): upper bound on how many times the tree may be split (**+** or **Auto**), so **k ≤ max splits + 1**.
  - **Min cluster size** (default **2**): any leaf with fewer than this many sequences is **noise** (−1 / `_cl-na`) and is excluded from the reported k.
  - **Remove long-branch tips as noise** (default **on**): when checked, at **open**, **Reset**, and the start of **Auto**, tips whose **incoming** branch length (`node.len`) is **greater than** the **Target max depth** field (linked to **≈ DPI** via **etd = 2 × dsr × DPI**) are marked **noise** and **excluded from splitting**. The remaining tips are re-rooted as one cluster; their tree-based DPI is recomputed (typically much lower if long tips dominated mean pairwise distance). Uncheck to cluster all tips.
  - **Target max depth** (default = **etd** from alignment load) with a linked **≈ DPI** field (rounded nearest day; either field updates the other via **etd = 2 × dsr × DPI**). **Auto** keeps splitting while some cluster’s **tree-based DPI** (mean pairwise path distance / `2×dsr`) is **greater than** that ≈ DPI target, and stops when all cluster DPIs ≤ target (or max splits). Bipartitions themselves still minimise **max(child LCA depths)** — Auto’s stop rule is DPI, not max depth. Editing target / DPI / min size / max splits does **not** clear the bubble tree; use **Reset** or **Auto**. Undersized children become **noise**. Auto may split pairs into noise-only leaves when needed.
- **Algorithm:**
  1. On open: download a CSV of the tip×tip **tree path-distance** matrix for all sample tips (symmetric, zero diagonal; filename `<alignment>_kdpis_tree_distances.csv`); optionally remove long-branch tips as pre-noise; build tree geometry for the remaining tips (splits + DPI).
  2. Start with one cluster (remaining sample tips; REF / SubType / PDB excluded). Leaf **max depth** = LCA→farthest-tip height; leaf **dpi** = mean pairwise tree path distance / (2 × dsr).
  3. On **+**: among leaves with size ≥ min size + 1 that are not marked unsplittable, split the one with the **largest max depth**. On **Auto**: among leaves with n ≥ 2 and **dpi** > target DPI, split by the same max-depth preference (noise-only children allowed).
  4. **Bipartition search:** consider every **internal phylogeny node** whose descendant tips form a proper nonempty subset of the cluster. Left = those tips; right = remainder. Score = **max(depth₁, depth₂)**. Prefer cuts where **both** sides meet min size, then one side, then (Auto only) neither. Choose the best score. If no valid bipartition exists, or the best score is **worse than** the parent’s max depth, mark the leaf **unsplittable** (still a normal cluster) and try the next candidate. If every remaining candidate is unsplittable, stop and **warn**.
  5. On **−**: reverse the most recent split. On **Reset**: clear splits, re-apply long-branch noise (if checked), return to one cluster of the remaining tips.
  6. Non-noise leaves are numbered **1…k by phylogeny tip order**: clusters nearer the **top** of the tree (smaller row index after ladderize) get **lower** IDs; those toward the **bottom** get **higher** IDs. Undersized leaves and pre-noise long-branch tips are **noise**.
- **Display:** Bubble **tree** of the split hierarchy: edges show parent→child splits; **leaf** bubbles use cluster colours (noise grey); internal nodes are grey. Radius ∝ **log(n)**; label **`n, dpi`** (tree mean-pairwise DPI).
- **Accept:** Writes `_cl-<id>` (or `_cl-na` for noise), stores **tree-based** cluster DPI / counts for the legend as **`id (n, dpi)`**, and **`Noise (n)`** when noise is present. **Cancel:** Clears preview clustering without keeping tags.

### 4.4 UMAP
- **Function:** Reduce pairwise leaf distances to 2D with **UMAP**, then cluster in 2D with DBSCAN (same as MDS after projection).
- **Algorithm:** UMAP (via `umap-js`) is run on the distance matrix (or a derived affinity matrix). Resulting 2D coordinates are then passed to the same DBSCAN + renumbering + Calinski–Harabasz pipeline as MDS.
- **Parameters:** **nNeighbors**, **Spread**, **Min Distance** (UMAP), plus **Radius (eps)** and **Min Neighbors** (DBSCAN). Eps max = half the larger projection axis range; step = max/50; min = step to avoid CH infinity; **Auto** on eps evaluates every slider tick and maximizes the selected index. **Estimate DPI** is available in the same dialog (see MDS).

### 4.5 MDS (Classical Multidimensional Scaling)
- **Function:** Reduce **pairwise leaf distances** (from the tree) to 2D, then cluster points in 2D with **DBSCAN**.
- **Projection:** Pairwise distances are taken from the tree (path length between leaves). **Classical MDS:** D² is double-centered to form **B** = −0.5 · H · D² · H, where **H** = I − (1/n)·1·1ᵀ (identity minus n⁻¹ times the matrix of ones). Top two eigenvalues and eigenvectors of B are computed (deterministic power iteration with fixed starting vectors and canonical axis signs); coordinates are the eigenvectors scaled by √λᵢ. Re-opening the MDS dialog on the same tree yields the same 2D layout.
- **Clustering:** DBSCAN on the 2D points (see §4.5). The subtype point is excluded from DBSCAN expansion and effectively not part of the density-based clustering (it is treated as an always-isolated reference). Clusters are renumbered by distance from the subtype reference (and, if a founder exists, may be used for ordering). Calinski–Harabasz is computed from the **tree** using the same cluster assignment (so all methods are comparable).
- **Parameters:** **Radius (eps)** and **Min Neighbors** for DBSCAN. Eps max = half the larger MDS axis range; step = max/50; min = step (avoids degenerate CH). **Auto** on eps: evaluates every slider tick, maximizes the selected index, updates slider and plot.
- **Estimate DPI:** Computes days-post-infection for each current (non-noise) cluster using the **mean** pairwise p-distance, then the same **JC69** correction and **dsr** as Group: **d = −(3/4) ln(1 − 4·api/3)**, **dpi = d / (2 × dsr)**. Shows results as **`clusterId (n, dpi)`** (singletons show `(n, —)`).
- **Apply:** Applies the clustering to the tree (writes cluster tags to names and tree), computes cluster DPI estimates with the current **dsr** (mean pairwise + JC69), and updates the legend as **`clusterId (n, dpi)`**. **Close:** Dismisses the overlay without applying.

### 4.6 DBSCAN (used in MDS/UMAP)
- **Algorithm:** Standard DBSCAN. Points within **eps** (Euclidean in 2D) are neighbors. If a point has ≥ **minPts** neighbors, it and all density-reachable points form a cluster. Otherwise it is noise (-1).

### 4.7. Calinski–Harabasz index (tree-based)

Used to score and optimize clusterings (tree-clade, tree-cut, MDS/UMAP). All use the **same** index on the **tree**.

- **Definitions:**  
  - **k** = number of clusters, **n** = number of leaves in those clusters.  
  - **Between-cluster (B):** For each cluster, take the MRCA of its leaves; *d*<sub>MRCA</sub> = branch distance from MRCA to root. Then B = Σ<sub>c</sub> n<sub>c</sub> · *d*<sub>MRCA</sub>² (for k ≥ 2; for k = 1 a synthetic B is used: (n/2)·(halfMax)² where halfMax = half the longest branch in the tree).  
  - **Within-cluster (W):** For each leaf, *d*<sub>leaf</sub> = distance to root; W = Σ over leaves of (*d*<sub>leaf</sub> − *d*<sub>MRCA</sub>)² for that leaf’s cluster.

- **Formula:**  
  **CH** = [ (B/(k−1)) / (W/(n−k)) ] × 1/(2<sup>k</sup>)  
  (for k = 1 the numerator is B and denominator W/(n−1); then the same 1/2<sup>k</sup> factor). The 2<sup>k</sup> term biases against large k.

- **Special cases:** If W ≤ 0, CH is returned as ∞ (best). Noise (−1) is excluded from clustering; remaining cluster IDs are used in the CH computation.

### 4.8. Ball–Hall-Adapted index (tree-based)

This index is displayed alongside the Calinski–Harabasz value in the tree-based hierarchical clustering dialogs.

- **Ball–Hall dispersion for the current clustering (k clusters):**  
  For each cluster c, let MRCA(c) be the most recent common ancestor of that cluster’s leaves, and let `dist(x, y)` be the branch-length distance between nodes. Each leaf contributes the squared distance from the leaf to its cluster’s MRCA:
  - `dist(leaf, MRCA(c)) = dist(leaf, root) − dist(MRCA(c), root)`
  - **BH(k) = Σ<sub>c</sub> Σ<sub>leaf in c</sub> dist(leaf, MRCA(c))²**

- **Ball–Hall dispersion for k=1 (whole tree):**  
  When there is only one cluster, MRCA(c) is the tree root, so:
  - **BH(1) = Σ<sub>leaf</sub> dist(leaf, root)²**

- **Ball–Hall-Adapted value (requested adaptation):**  
  The adapted index compares dispersion for the full tree to dispersion for the current clustering, then applies the same **2<sup>k</sup>** penalty used for Calinski–Harabasz:
  - **BH<sub>adapted</sub> = [ BH(1) / BH(k) ] / 2<sup>k</sup>**

If BH(k) = 0, the adapted index is returned as ∞.

---

## 5. Epitopes

### Load Epitopes
- **Function:** Load epitope definitions from a **CSV** file.
- **CSV format:** One row per epitope. First column = epitope **name**. Remaining columns = regions: either `start:end` or a single position (same start and end). Example:

  `VRC01,197:198,230,276,278:282,365:371,427:428,430,455:463,465,469,471:474`  
  `CAP256,156:163,166:167,169:170,178:179,181:184`
- **Result:** Epitopes are stored in `state.epitopes`. If the alignment was in NT mode, it is converted to AA (and mode switched to AA). Select Epitope and Export Epitopes are enabled.
- **Merge behavior:** Loading new epitopes **merges** with existing ones: same name overwrites, other names are kept.

### Select Epitope
- **Function:** Choose which epitope is active for display (e.g. alignment columns and logo). Option “None” shows all columns.
- **Usage:** Opens a dropdown of loaded epitopes; selection restricts the visible alignment to that epitope’s regions.

### Show Logo
- **Function:** Generate a **sequence logo** for the currently selected epitope (and selected sequences).
- **Requirement:** An epitope must be selected. Logo shows conservation per position in the epitope regions.

### Export Epitopes
- **Function:** Download the current set of epitopes as **epitopes.csv**.
- **Format:** Same CSV as load (name, then region columns as `start:end` or single positions). Only epitopes with at least one region are exported.
- **State:** Disabled when there are no epitopes.

### New epitope (3D)
- **Function:** Define an epitope from the loaded 3D structure by **distance**: residues on an “intrinsic” chain whose minimum distance to an “extrinsic” chain is below a cutoff.
- **Parameters:** Intrinsic chain, extrinsic chain, **distance cutoff (Å)**. Name is auto-generated as `<pdbCode>_<cutoff>A`.
- **Result:** Residues are mapped to alignment columns (via reference); epitope is added to `state.epitopes` and selected. Select Epitope and Export Epitopes are enabled if they were not already.

---

## 6. PDB structure

### Load file / Fetch
- **Function:** Load a PDB/CIF from file or fetch by **PDB ID** from RCSB. Chains are added as sequences (e.g. `xxxx_Chain_A [PDB_A]`) and aligned to the reference in AA mode.
- **Usage:** “Choose file” or enter a 4-character ID and click Fetch.

### View 3D
- **Function:** Open the 3D viewer (3Dmol.js). Chains can be toggled and styled. Epitope residues (from loaded or 3D-defined epitopes) are colored on intrinsic chains; rest can be gray.
- **Chain categories:** The 3D viewer categorizes chains into three groups:
  - **Intrinsic (≥50%)**: Chains present in the alignment panel that align at **≥50%** coverage/identity to the alignment reference.
  - **Extrinsic (<50%)**: Chains present in the alignment panel that align at **<50%** coverage/identity.
  - **Superficial**: Chains **not** present in the alignment panel.
- **Epitope residue coloring** is applied to intrinsic chains; clicking on an active epitope residue initiates a **partial logo** dialogue; you can show/hide chains and set styles separately per category.
- **New epitope:** Define a distance-based epitope from the 3D structure (see §5).
- **Other buttons:** Full screen, Export PNG, Close.

---

## 7. Other controls

### fasta (Download)
- **Function:** Download the **current view** (NT or AA) as FASTA. The download uses the **uploaded alignment filename**. If clustering is active, **`_clus`** is inserted before the extension (e.g. `myfile.fasta` → `myfile_clus.fasta`).

### Clear Alignment
- **Function:** Remove all sequences and reset state (groups, clusters, tree, founder designation, etc.). Returns to initial empty state.

### Help (?)
- **Function:** Toggle the help overlay with short descriptions of each control. For full documentation, open the published **aliViz user guide** (link in the help overlay).

---

## 8. Color legend

- **Groups:** Lists group labels and colors for sequence **names** (from Group). When DPI estimates are available, labels include sequence count and JC69-based dpi as **`label (n, dpi)`**, e.g. **`2000 (12, 56)`**.
- **Clusters:** Lists cluster IDs (numeric labels) and colors for tree tips (from any clustering method). After **Accept**/**Apply**, labels include **`clusterId (n, dpi)`** when available, e.g. **`1 (8, 56)`**, using the same JC69 + **dsr** formula as groups. Noise appears when applicable. **Red is not used** in the cluster palette (reserved for the non-functional leaf oval on trees).
- **Note:** The floating legend does **not** include an alignment-length or phylogenetic scale bar; phylogenetic scale bars appear only on **exported SVG** tree figures (see §3).
- Cluster colors are removed when clustering is cleared (e.g. **None** or **Cancel**). Group colors are preserved after Cancel by rekeying group membership to the stripped (no `_cl-*`) names.

---

## 9. Summary of algorithms and formulas

| Feature        | Algorithm / formula |
|----------------|---------------------|
| Load sanitize  | A–Z, `-`, `*` kept; other chars → `X`; alert user. |
| Functionality  | AA: ungapped M start, `*` end, no internal `*`; NT: translate frame 1 then same test; tag `functional`; report `s,f,nf,dpi,dsr,etd` above filename (all seqs; DPI from filename `dpi±N` or 14; etd=2×dsr×DPI) and `s,f,nf` in SVG title (tree leaves only); red oval around non-functional tree tips. |
| Has subtype    | If off, no subtype index; toggling clears group/tree/cluster state. |
| Group          | Split name by delimiter; group = field value; unique IDs; dpi via mean pairwise p-distance + JC69: d=−(3/4)ln(1−4·api/3), dpi=d/(2×dsr) (dsr from Group dialog or on-load DPI prompt, default 7.9e-5 for the session); legend `label (n, dpi)`. |
| Sort           | REF, SubType (if any) fixed; others by group ID then name. |
| Prune          | Remove sequences in selected groups; renumber group IDs. |
| Founder consensus | Majority per column over selected group; `consensus_of_*` name. |
| Founder medoid | Minimize sum of padded Hamming distances to other group sequences; keep original name; `[Founder]` label. |
| NJ tree        | PhyloTools from alignment (pairwise distances → NJ); optional auto root on founder + ladderize by depth. |
| FastTree       | bioWASM FastTree on alignment; optional auto root on founder + ladderize by depth. |
| Load Newick    | Parse Newick; validate leaf names vs alignment; optional auto root on founder + ladderize by depth (Infer dialog). |
| Reroot         | Find founder leaf; reroot on edge to that leaf. |
| Ladderize      | By weight: sort children by leaf count. By depth: sort by max root-to-leaf depth in subtree. |
| Histogram      | Root-to-leaf distance = sum of branch lengths; bin and plot. |
| SVG scale      | Auto: fit to 520 px (linear) / R=260 (circular, ½ px per unit). Fixed: branch length per cm; overflow check. DPI: days default from filename (`dpi±N`, max if several) else 14; MPD seeded from load-time dsr; expected root→tip depth (DPI × MPD) → 32.5 px = 1/16 mark (linear) / 16.25 px (circular); no overflow abort (truncate at right border, half diamond + red oval if non-functional on clipped tips); dotted expected-depth line. |
| Tree-clade     | DFS; new cluster when edge length > threshold; min-leaves → noise. |
| Tree-cut       | Three depth cutoffs; BFS assign clusters; min-leaves → noise. |
| k-DPIs         | Divisive tree: split at internal node minimising max(depth₁,depth₂) (LCA height); optional pre-noise if incoming branch &gt; target etd (default on); failed improve → unsplittable (kept as cluster), try next; if all unsplittable → warn/stop; +/−; Max splits (default 9); Reset → k=1; Auto until all cluster DPIs ≤ target (tree mean pairwise / 2dsr) or max splits; min size → noise; tree mean-pairwise DPI labels; CSV of tip×tip tree distances on open; bubble tree. |
| Cluster None   | Clear `leafClusters`; strip `_cl-*`; keep groups. |
| MDS            | B = −0.5·H·D²·H; eigendecomposition; coords = eigenvectors × √(eigenvalues). |
| UMAP           | External UMAP on distances → 2D. |
| DBSCAN         | eps-neighborhood; minPts; density-connected components; subtype isolated. |
| Calinski–Harabasz | B/(k−1) and W/(n−k) on tree (MRCA/root distances); divide by 2<sup>k</sup>. |
| Epitope CSV    | Rows: name, region1, region2, … (e.g. `start:end` or single position). |
| Distance epitope | Min distance intrinsic→extrinsic &lt; cutoff (Å); map residues to alignment. |

---

*End of aliViz User Guide.*
