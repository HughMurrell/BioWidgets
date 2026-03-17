# aliViz User Guide

aliViz is a bioinformatics alignment and phylogeny viewer. It supports loading alignments, inferring and manipulating trees, grouping sequences, clustering (tree-based and projection-based), epitope definition and logos, and 3D structure viewing with epitope coloring.

---

## 1. Loading and view controls

### Choose File
- **Function:** Load an alignment from a file (FASTA or FASTQ). The first sequence is treated as **Reference (REF)**, the second as **SubType**. Remaining sequences are aligned and displayed.
- **Usage:** Click the file input or “Choose File” to select an alignment file.

### View Mode: NT / AA
- **Function:** Switch between **Nucleotide (NT)** and **Amino acid (AA)** view.
- **Algorithm:** In AA mode, nucleotides are translated using the standard genetic code (frame can be set if frame selector is enabled). Gaps and invalid codons produce gap characters or `X`.
- **Modes:** **NT** (default), **AA**.

### AA palette
- **Function:** Choose the amino acid color scheme.
- **Options:** *Alignment (IUPAC)* (default) or *Biochemical* (e.g. by property).

### Highlighter
- **Function:** Highlight differences in the alignment with respect to a chosen reference.
- **Options:** **Off** (default), **Founder** (consensus of a group), **SubType**, **Reference**. Mismatches against the selected sequence are highlighted.

---

## 2. Grouping and sequences

### Group
- **Function:** Assign sequences to groups using a delimiter and a field in the sequence name.
- **Algorithm:** You specify a **delimiter** (e.g. `_`) and a **field number** (1-based). Each sequence name is split by the delimiter; the value at that field becomes the group label. REF, SubType, and PDB chain sequences are excluded from grouping. Groups are assigned unique IDs and used for coloring and legend.
- **Result:** `state.sequenceGroups` (name → group ID) and `state.groupNames` (group ID → label) are set; sequence names are colored by group in the name panel.

### Sort
- **Function:** Reorder sequences by current grouping (or by name if no grouping), keeping REF first and SubType second.
- **Algorithm:** REF and SubType stay at the top; other sequences are sorted by group ID (from `state.sequenceGroups`), then by name within each group. PDB chains and founder are handled in the sort order.

### Prune
- **Function:** Remove sequences that belong to selected groups.
- **Usage:** Open the Prune overlay, select groups to **remove**, then Apply. Removed sequences are deleted from the alignment. Group IDs are renumbered to 0, 1, 2, … for the remaining groups.

### Add Founder
- **Function:** Add a **founder** (consensus) sequence for a chosen group.
- **Algorithm:** You select a group. A consensus sequence is built from that group (e.g. majority character per column, gaps where no majority). It is inserted into the alignment and named e.g. `consensus_of_<groupLabel>`. The founder can be used as the Highlighter reference and as a reroot target.

---

## 3. Tree inference and manipulation

### Infer
- **Function:** Infer a phylogeny from the current alignment (excluding REF, SubType, and PDB chains).
- **Methods:**
  - **FastTree** (default): Uses FastTree (via bioWASM) to infer a maximum-likelihood style tree from the alignment.
  - **Neighbor Joining (NJ):** Uses PhyloTools to compute pairwise distances and build an NJ tree from the alignment.
- **Note:** Before inference, any existing clustering is cleared and cluster suffixes (`_cl-*`) are removed from sequence and tree node names.

### Reroot
- **Function:** Reroot the tree on the **Founder** sequence.
- **Targets:** **Founder** (default). *(The current app requires a founder sequence for rerooting.)*
- **Algorithm:** The founder leaf is located and the tree is rerooted at the edge leading to it so that the founder is the outgroup. Branch lengths and topology are preserved; only the root position changes.

### Ladderize
- **Function:** Reorder the tree and the sequence list so that the tree is ladderized and sequences match leaf order.
- **Modes:**
  - **By Depth** (default): At each node, compute the maximum root-to-leaf distance in each child’s subtree (using branch lengths). Sort children by that depth (ascending). Shallower subtrees appear first.
  - **By Weight:** At each node, sort children by number of leaves (ascending). Lighter subtrees appear first.
- **Result:** Tree drawing order and `state.viewSequences` are updated to follow the new leaf order. Histogram, Cluster, and tree export are enabled after ladderizing.

### Histogram
- **Function:** Show a histogram of **root-to-leaf distances** (sum of branch lengths from root to each leaf).
- **Algorithm:** For each leaf, the path from root to leaf is traversed and branch lengths (`node.len`) are summed. Distances are binned (number of bins between 10 and 30, about √n). Bars are drawn for each bin count.
- **Requirement:** Tree must be loaded and ladderized.

### tree (Export tree)
- **Function:** Export the current tree in **Newick** format.
- **Usage:** Click “tree” to download a `.txt` or Newick file containing the current tree with names and branch lengths (if present).

### Clear Tree
- **Function:** Remove the current tree. Clustering state is cleared. Histogram, Cluster, and related buttons are disabled.

---

## 4. Clustering (Cluster button)

Clicking **Cluster** opens a method selector, then the corresponding clustering interface.
Each method allows the user to cluster sequences. Special sequences such as Reference, Subtype or any PDB chains are excluded from the clustering.
- **Methods (selector):** **Hierarchical (tree-clade clustering)** (default), **Hierarchical (tree-cut clustering)**, **UMAP**, **MDS**.

### 4.1 Hierarchical (tree-clade clustering)
- **Function:** Define clusters by **clades** on the tree: a new clade (and thus a new cluster) starts when the **incoming branch length** to a node exceeds a threshold.
- **Parameters:**
  - **Branch length threshold:** Minimum branch length that starts a new clade (slider from step to max). Leaves are assigned to the clade that contains them.
  - **Min leaves per cluster:** Clusters with fewer than this many leaves are relabelled as noise (`_cl-na`).
- **Algorithm:** DFS from root. When traversing an edge longer than the threshold, increment cluster ID. Assign each leaf to the current cluster. Then apply min-leaves filter, mark small clusters as noise, and renumber clusters 1, 2, 3, …
- **Auto:** Searches the threshold (from step to max) to **maximize the Calinski–Harabasz index** (see §6) for the current min-leaves. Sets the slider to the best value.
- **Accept:** Applies the clustering: adds `_cl-<id>` or `_cl-na` to sequence and tree node names, updates the legend and tree colors.
- **Cancel:** Removes all clustering: clears `state.leafClusters`, strips `_cl-*` from names and tree nodes, updates legend and redraws. **Group colors** on sequence names are preserved (group mapping is rekeyed to the stripped names).

### 4.2 Hierarchical (tree-cut clustering)
- **Function:** Cut the tree with **three depth lines** (root-to-node distance). Each region between/above lines defines clusters.
- **Parameters:**
  - **Level 1, 2, 3 Depth:** Sliders (0 to max root-to-leaf distance). Order is enforced: Level 1 ≤ Level 2 ≤ Level 3.
  - **Min leaves per cluster:** As in tree-clade; clusters below this size become noise.
- **Algorithm:** Leaves with depth &lt; Level 1 = cluster 1. BFS from root: when a node’s depth crosses Level 1 (or 2 or 3), all leaves in its subtree with depth ≥ that level are assigned the next cluster ID. Then min-leaves → noise, renumber clusters.
- **Auto:** **Full grid search** over all triples (d1, d2, d3) with d1 ≤ d2 ≤ d3 (sampled from slider range). For each triple, computes the cluster map and the Calinski–Harabasz index; chooses the triple that maximizes it and sets the three sliders.
- **Accept / Cancel:** Same idea as tree-clade: Accept writes cluster tags and updates legend/tree; Cancel clears clustering and strips `_cl-*` while keeping group colors.

### 4.3 UMAP
- **Function:** Reduce pairwise leaf distances to 2D with **UMAP**, then cluster in 2D with DBSCAN (same as MDS after projection).
- **Algorithm:** UMAP (via `umap-js`) is run on the distance matrix (or a derived affinity matrix). Resulting 2D coordinates are then passed to the same DBSCAN + renumbering + Calinski–Harabasz pipeline as MDS.
- **Parameters:** **nNeighbors**, **Spread**, **Min Distance** (UMAP), plus **Radius (eps)** and **Min Neighbors** (DBSCAN). Eps min = step to avoid CH infinity; **Auto** on eps optimizes CH over the eps range.

### 4.4 MDS (Classical Multidimensional Scaling)
- **Function:** Reduce **pairwise leaf distances** (from the tree) to 2D, then cluster points in 2D with **DBSCAN**.
- **Projection:** Pairwise distances are taken from the tree (path length between leaves). **Classical MDS:** D² is double-centered to form **B** = −0.5 · H · D² · H, where **H** = I − (1/n)·1·1ᵀ (identity minus n⁻¹ times the matrix of ones). Top two eigenvalues and eigenvectors of B are computed; coordinates are the eigenvectors scaled by √λᵢ.
- **Clustering:** DBSCAN on the 2D points (see §4.5). The subtype point is excluded from DBSCAN expansion and effectively not part of the density-based clustering (it is treated as an always-isolated reference). Clusters are renumbered by distance from the subtype reference (and, if a founder exists, may be used for ordering). Calinski–Harabasz is computed from the **tree** using the same cluster assignment (so all methods are comparable).
- **Parameters:** **Radius (eps)** and **Min Neighbors** for DBSCAN. Eps slider range is step to max (no 0) to avoid degenerate CH. **Auto** on eps: varies eps from step to max, maximizes CH, updates slider and plot.
- **Apply:** Applies the clustering to the tree (writes cluster tags to names and tree, updates legend). **Close:** Dismisses the overlay without applying.

### 4.5 DBSCAN (used in MDS/UMAP)
- **Algorithm:** Standard DBSCAN. Points within **eps** (Euclidean in 2D) are neighbors. If a point has ≥ **minPts** neighbors, it and all density-reachable points form a cluster. Otherwise it is noise (-1).

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

## 6. Calinski–Harabasz index (tree-based)

Used to score and optimize clusterings (tree-clade, tree-cut, MDS/UMAP). All use the **same** index on the **tree**.

- **Definitions:**  
  - **k** = number of clusters, **n** = number of leaves in those clusters.  
  - **Between-cluster (B):** For each cluster, take the MRCA of its leaves; *d*<sub>MRCA</sub> = branch distance from MRCA to root. Then B = Σ<sub>c</sub> n<sub>c</sub> · *d*<sub>MRCA</sub>² (for k ≥ 2; for k = 1 a synthetic B is used: (n/2)·(halfMax)² where halfMax = half the longest branch in the tree).  
  - **Within-cluster (W):** For each leaf, *d*<sub>leaf</sub> = distance to root; W = Σ over leaves of (*d*<sub>leaf</sub> − *d*<sub>MRCA</sub>)² for that leaf’s cluster.

- **Formula:**  
  **CH** = [ (B/(k−1)) / (W/(n−k)) ] × 1/(2<sup>k</sup>)  
  (for k = 1 the numerator is B and denominator W/(n−1); then the same 1/2<sup>k</sup> factor). The 2<sup>k</sup> term biases against large k.

- **Special cases:** If W ≤ 0, CH is returned as ∞ (best). Cluster 0 (subtype) and noise (−1) are included in the tree-based computation where applicable.
- **Special cases:** If W ≤ 0, CH is returned as ∞ (best). Noise (−1) is excluded from clustering; remaining cluster IDs are used in the CH computation.

---

## 7. PDB structure

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

## 8. Other controls

### fasta (Download)
- **Function:** Download the **current view** (NT or AA) as FASTA. Filename includes mode, e.g. `alignment_AA.fasta`.

### Clear Alignment
- **Function:** Remove all sequences and reset state (groups, clusters, tree, etc.). Returns to initial empty state.

### Help (?)
- **Function:** Toggle the help overlay with short descriptions of each control.

---

## 9. Color legend

- **Groups:** Lists group labels and colors for sequence **names** (from Group).
- **Clusters:** Lists cluster IDs for cluster tags and colors for tree tips (from any clustering method).  
Cluster colors are removed when clustering is cleared (e.g. Cancel). Group colors are preserved after Cancel by rekeying group membership to the stripped (no `_cl-*`) names.

---

## 10. Summary of algorithms and formulas

| Feature        | Algorithm / formula |
|----------------|---------------------|
| Group          | Split name by delimiter; group = field value; unique IDs. |
| Sort           | REF, SubType fixed; others by group ID then name. |
| Prune          | Remove sequences in selected groups; renumber group IDs. |
| Founder        | Consensus (e.g. majority) per column over selected group. |
| NJ tree        | PhyloTools from alignment (pairwise distances → NJ). |
| FastTree       | bioWASM FastTree on alignment. |
| Reroot         | Find target leaf; reroot on edge to that leaf. |
| Ladderize      | By weight: sort children by leaf count. By depth: sort by max root-to-leaf depth in subtree. |
| Histogram      | Root-to-leaf distance = sum of branch lengths; bin and plot. |
| Tree-clade     | DFS; new cluster when edge length > threshold; min-leaves → noise. |
| Tree-cut       | Three depth cutoffs; BFS assign clusters; min-leaves → noise. |
| MDS            | B = −0.5·H·D²·H; eigendecomposition; coords = eigenvectors × √(eigenvalues). |
| UMAP           | External UMAP on distances → 2D. |
| DBSCAN         | eps-neighborhood; minPts; density-connected components; subtype isolated. |
| Calinski–Harabasz | B/(k−1) and W/(n−k) on tree (MRCA/root distances); divide by 2<sup>k</sup>. |
| Epitope CSV    | Rows: name, region1, region2, … (e.g. `start:end` or single position). |
| Distance epitope | Min distance intrinsic→extrinsic &lt; cutoff (Å); map residues to alignment. |

---

*End of aliViz User Guide.*
