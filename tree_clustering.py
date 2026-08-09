#!/usr/bin/env python3
"""
Batch k-DPIs clustering + SVG export (aliViz-compatible).

For each Newick tree in an input directory:
  1. Parse DPI from the filename (dpi±N; largest N; default 14)
  2. Group tip names by delimiter/field (defaults match aliViz: '_' field 3)
  3. Pick the smallest group; choose a founder as the tree-path medoid
  4. Reroot on the founder; ladderize by depth (shallow subtrees first)
  5. Run k-DPIs Auto clustering (same bipartition / DPI rules as aliViz)
  6. Write a linear SVG (default: AUTO scale with DPI overlay when dpi±N is in
     the filename; legacy DPI 1/16 scale via --scale dpi)
  7. After the directory is processed, write housekeeping.csv (CAPid, Groups, Clusters)
     and housekeeping.svg (pie chart by cluster-count bins: 1, 2, 3–4, 5–8, …)

Usage:
  python tree_clustering.py INPUT_DIR OUTPUT_DIR
  python tree_clustering.py INPUT_DIR OUTPUT_DIR --delimiter _ --field 3
  python tree_clustering.py INPUT_DIR OUTPUT_DIR --scale dpi   # legacy DPI scale
"""

from __future__ import annotations

import argparse
import csv
import math
import re
import sys
from collections import defaultdict
from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterable, Optional


# --- Constants matching aliViz.html ---
DEFAULT_DSR = 0.00005
DEFAULT_DPI_DAYS = 14
DEFAULT_DELIMITER = "_"
DEFAULT_FIELD = 3  # 1-based, same as aliViz Group dialog
DEFAULT_MIN_SIZE = 1
DEFAULT_MAX_SPLITS = 63
DEFAULT_SCALE_MODE = "auto"  # "auto" (aliViz AUTO) or "dpi" (legacy 1/16 mark)
BRANCH_W = 520.0
DPI_MARK_PX = BRANCH_W / 16.0  # 32.5 — used only by legacy --scale dpi
ROW_H = 20.0
NOISE_COLOR = "#9ca3af"
COLOR_PALETTE = [
    "#ff00ff",
    "#3b82f6",
    "#10b981",
    "#f59e0b",
    "#8b5cf6",
    "#06b6d4",
    "#84cc16",
    "#f97316",
    "#6366f1",
    "#14b8a6",
    "#64748b",
    "#fbbf24",
    "#a855f7",
    "#22c55e",
]

TREE_EXTENSIONS = {".nwk", ".newick", ".tree", ".tre", ".treefile", ".txt"}
DPI_RE = re.compile(r"dpi[+-](\d+)", re.IGNORECASE)
CAP_RE = re.compile(r"(CAP\d+)", re.IGNORECASE)
CLUSTER_SUFFIX_RE = re.compile(r"_(?:cl-(?:\d+|na|n))$", re.IGNORECASE)


# ---------------------------------------------------------------------------
# Tree structure
# ---------------------------------------------------------------------------

@dataclass(eq=False)
class Node:
    name: Optional[str] = None
    length: float = 0.0  # incoming branch length
    children: list["Node"] = field(default_factory=list)
    parent: Optional["Node"] = None
    # layout
    x_depth: float = 0.0
    y_row: float = 0.0
    subtree_max_depth: float = 0.0

    def __hash__(self) -> int:
        return id(self)

    @property
    def is_leaf(self) -> bool:
        return not self.children

    def leaves(self) -> list["Node"]:
        if self.is_leaf:
            return [self]
        out: list[Node] = []
        for c in self.children:
            out.extend(c.leaves())
        return out


@dataclass
class Tree:
    root: Node

    def tips(self) -> list[Node]:
        return self.root.leaves()

    def tip_names(self) -> list[str]:
        return [t.name for t in self.tips() if t.name]


# ---------------------------------------------------------------------------
# Newick I/O
# ---------------------------------------------------------------------------

def _parse_newick_subtree(s: str, i: int) -> tuple[Node, int]:
    node = Node()
    if i < len(s) and s[i] == "(":
        i += 1
        while True:
            child, i = _parse_newick_subtree(s, i)
            node.children.append(child)
            child.parent = node
            if i >= len(s):
                raise ValueError("Unterminated Newick subtree")
            if s[i] == ",":
                i += 1
                continue
            if s[i] == ")":
                i += 1
                break
            raise ValueError(f"Unexpected '{s[i]}' in Newick at {i}")
    # name
    start = i
    while i < len(s) and s[i] not in "(),:;[":
        i += 1
    name = s[start:i].strip()
    if name:
        node.name = name
    # length
    if i < len(s) and s[i] == ":":
        i += 1
        start = i
        while i < len(s) and s[i] in "0123456789.eE+-":
            i += 1
        try:
            node.length = float(s[start:i])
        except ValueError:
            node.length = 0.0
    return node, i


def parse_newick(text: str) -> Tree:
    s = "".join(text.split())
    if not s:
        raise ValueError("Empty Newick")
    if s[-1] == ";":
        s = s[:-1]
    root, i = _parse_newick_subtree(s, 0)
    if i != len(s):
        # trailing junk ignored if only whitespace already stripped
        pass
    return Tree(root=root)


def to_newick(node: Node) -> str:
    if node.is_leaf:
        name = node.name or ""
        return f"{name}:{node.length:.8g}" if node.length else name
    inner = ",".join(to_newick(c) for c in node.children)
    name = node.name or ""
    if node.length:
        return f"({inner}){name}:{node.length:.8g}"
    return f"({inner}){name}"


# ---------------------------------------------------------------------------
# Tree geometry helpers
# ---------------------------------------------------------------------------

def set_parents(node: Node, parent: Optional[Node] = None) -> None:
    node.parent = parent
    for c in node.children:
        set_parents(c, node)


def node_to_root_distance(node: Node) -> dict[Node, float]:
    """Distance from tree root along branches to each node (includes node.length)."""
    dist: dict[Node, float] = {}

    def walk(n: Node, d: float) -> None:
        cur = d + (n.length or 0.0)
        dist[n] = cur
        for c in n.children:
            walk(c, cur)

    # Caller's root should already be the phylogeny root; its length is usually 0.
    walk(node, 0.0)
    return dist


def ancestors(node: Node) -> list[Node]:
    out = []
    cur: Optional[Node] = node
    while cur is not None:
        out.append(cur)
        cur = cur.parent
    return out


def lca(a: Node, b: Node) -> Node:
    aset = set(ancestors(a))
    cur: Optional[Node] = b
    while cur is not None:
        if cur in aset:
            return cur
        cur = cur.parent
    return a  # fallback


def path_distance(a: Node, b: Node, root_dist: dict[Node, float]) -> float:
    if a is b:
        return 0.0
    anc = lca(a, b)
    return root_dist[a] + root_dist[b] - 2.0 * root_dist[anc]


def max_depth_of_tips(tips: list[Node], root_dist: dict[Node, float]) -> float:
    """LCA-subtree height: max path from LCA of tips to any tip."""
    if len(tips) < 2:
        return 0.0
    anc = tips[0]
    for t in tips[1:]:
        anc = lca(anc, t)
    base = root_dist[anc]
    return max((root_dist[t] - base) for t in tips)


def mean_pairwise_path(tips: list[Node], root_dist: dict[Node, float]) -> Optional[float]:
    n = len(tips)
    if n < 2:
        return None
    total = 0.0
    pairs = 0
    for i in range(n):
        for j in range(i + 1, n):
            total += path_distance(tips[i], tips[j], root_dist)
            pairs += 1
    return total / pairs if pairs else None


def dpi_days_from_mean(mean_d: Optional[float], dsr: float = DEFAULT_DSR) -> Optional[int]:
    if mean_d is None or not (dsr > 0):
        return None
    dpi = mean_d / (2.0 * dsr)
    if not math.isfinite(dpi) or dpi < 0:
        return None
    return int(round(dpi))


# ---------------------------------------------------------------------------
# Reroot & ladderize
# ---------------------------------------------------------------------------

def reroot_on_leaf(tree: Tree, leaf: Node) -> Tree:
    """Reroot so `leaf` is an outgroup child of a new root (edge to leaf is bisected)."""
    set_parents(tree.root)
    if leaf.parent is None:
        return tree

    path: list[Node] = []
    cur: Optional[Node] = leaf
    while cur is not None:
        path.append(cur)
        cur = cur.parent
    # path[0]=leaf, path[1]=parent, ..., path[-1]=old_root

    parent = path[1]
    parent.children = [c for c in parent.children if c is not leaf]
    leaf.parent = None
    orig_lens = [n.length or 0.0 for n in path]
    leaf_edge = orig_lens[0]

    # Invert the spine parent → … → old_root so `parent` roots the remaining tree.
    for i in range(1, len(path) - 1):
        node = path[i]
        up = path[i + 1]
        up.children = [c for c in up.children if c is not node]
        if up not in node.children:
            node.children.append(up)
        # Edge that was up→node (length on node) becomes node→up after flip.
        up.length = orig_lens[i]
        up.parent = node

    rest = parent
    new_root = Node(name=None, length=0.0)
    half = leaf_edge / 2.0
    leaf.length = half
    rest.length = half
    leaf.parent = new_root
    rest.parent = new_root
    new_root.children = [leaf, rest]
    set_parents(new_root)
    return Tree(root=new_root)


def ladderize_by_depth(tree: Tree) -> None:
    """Sort children by subtree max root→leaf depth ascending (aliViz depth mode)."""

    def compute(n: Node) -> float:
        if n.is_leaf:
            n.subtree_max_depth = n.length or 0.0
            return n.subtree_max_depth
        best = 0.0
        for c in n.children:
            d = compute(c)
            best = max(best, d)
        n.subtree_max_depth = (n.length or 0.0) + best
        return n.subtree_max_depth

    def sort_rec(n: Node) -> None:
        if n.is_leaf:
            return
        for c in n.children:
            sort_rec(c)
        n.children.sort(key=lambda c: c.subtree_max_depth)

    compute(tree.root)
    sort_rec(tree.root)
    set_parents(tree.root)


def layout_tree(tree: Tree) -> float:
    """Assign x_depth and y_row; return maxDepth."""
    set_parents(tree.root)

    def set_x(n: Node, depth: float) -> float:
        n.x_depth = depth + (n.length or 0.0)
        mx = n.x_depth
        for c in n.children:
            mx = max(mx, set_x(c, n.x_depth))
        return mx

    max_d = set_x(tree.root, 0.0)

    tips = tree.tips()
    for i, t in enumerate(tips):
        t.y_row = float(i)

    def set_y(n: Node) -> float:
        if n.is_leaf:
            return n.y_row
        ys = [set_y(c) for c in n.children]
        n.y_row = sum(ys) / len(ys) if ys else 0.0
        return n.y_row

    set_y(tree.root)
    return max_d


# ---------------------------------------------------------------------------
# Filename / grouping
# ---------------------------------------------------------------------------

def detect_dpi_from_filename(name: str) -> Optional[int]:
    """Largest dpi±N in filename, or None if absent."""
    matches = [int(m) for m in DPI_RE.findall(name)]
    return max(matches) if matches else None


def parse_dpi_from_filename(name: str) -> int:
    detected = detect_dpi_from_filename(name)
    return detected if detected is not None else DEFAULT_DPI_DAYS


def parse_capid_from_filename(name: str) -> str:
    """Return CAP token from filename (e.g. CAP008), or empty string if absent."""
    m = CAP_RE.search(name)
    return m.group(1).upper() if m else ""


def capid_sort_key(capid: str) -> tuple[int, str]:
    """Numeric CAP sort with missing ids last."""
    m = re.search(r"(\d+)", capid or "")
    if not m:
        return (10**12, capid or "")
    return (int(m.group(1)), capid)


def normalize_tip_name(name: str) -> str:
    return CLUSTER_SUFFIX_RE.sub("", name or "")


def group_tips(
    tip_names: Iterable[str],
    delimiter: str = DEFAULT_DELIMITER,
    field_1based: int = DEFAULT_FIELD,
) -> tuple[dict[str, int], dict[int, str], dict[int, list[str]]]:
    """Return name→groupId, groupId→label, groupId→member names. Excludes tips lacking the field."""
    labels_by_name: dict[str, str] = {}
    for raw in tip_names:
        nm = normalize_tip_name(raw)
        parts = nm.split(delimiter)
        idx = field_1based - 1
        if 0 <= idx < len(parts) and parts[idx] != "":
            labels_by_name[nm] = parts[idx]

    unique = sorted(set(labels_by_name.values()))
    label_to_id = {lab: i for i, lab in enumerate(unique)}
    name_to_gid = {nm: label_to_id[lab] for nm, lab in labels_by_name.items()}
    gid_to_label = {i: lab for lab, i in label_to_id.items()}
    members: dict[int, list[str]] = defaultdict(list)
    for nm, gid in name_to_gid.items():
        members[gid].append(nm)
    return name_to_gid, gid_to_label, dict(members)


def smallest_group_id(members: dict[int, list[str]]) -> Optional[int]:
    if not members:
        return None
    return min(members.keys(), key=lambda g: (len(members[g]), g))


def tree_medoid(tip_nodes: list[Node], root_dist: dict[Node, float]) -> Node:
    """Tip minimizing sum of tree path distances to other tips in the set."""
    if len(tip_nodes) == 1:
        return tip_nodes[0]
    best = tip_nodes[0]
    best_sum = float("inf")
    for i, a in enumerate(tip_nodes):
        s = 0.0
        for j, b in enumerate(tip_nodes):
            if i == j:
                continue
            s += path_distance(a, b, root_dist)
        if s < best_sum:
            best_sum = s
            best = a
    return best


# ---------------------------------------------------------------------------
# k-DPIs
# ---------------------------------------------------------------------------

@dataclass
class KDpisLeaf:
    tips: list[Node]
    dpi: Optional[int]
    depth: float
    unsplittable: bool = False
    left: Optional["KDpisLeaf"] = None
    right: Optional["KDpisLeaf"] = None
    cluster_id: int = 0
    is_noise: bool = False

    @property
    def n(self) -> int:
        return len(self.tips)

    @property
    def is_atomic(self) -> bool:
        return self.left is None and self.right is None


def make_kdpis_leaf(tips: list[Node], root_dist: dict[Node, float], dsr: float) -> KDpisLeaf:
    mean_d = mean_pairwise_path(tips, root_dist)
    return KDpisLeaf(
        tips=list(tips),
        dpi=dpi_days_from_mean(mean_d, dsr),
        depth=max_depth_of_tips(tips, root_dist),
    )


def collect_kdpis_leaves(node: KDpisLeaf, out: Optional[list[KDpisLeaf]] = None) -> list[KDpisLeaf]:
    acc = out if out is not None else []
    if node.is_atomic:
        acc.append(node)
        return acc
    assert node.left and node.right
    collect_kdpis_leaves(node.left, acc)
    collect_kdpis_leaves(node.right, acc)
    return acc


def _is_better_bip(cand: tuple[int, float], best: Optional[tuple[int, float]]) -> bool:
    if best is None:
        return True
    if cand[0] < best[0]:
        return True
    if cand[0] > best[0]:
        return False
    return cand[1] + 1e-9 < best[1]


def best_bipartition(
    tips: list[Node],
    tree_root: Node,
    root_dist: dict[Node, float],
    min_size: int,
) -> Optional[tuple[list[Node], list[Node], float]]:
    n = len(tips)
    if n < 2 * min_size:
        return None
    tip_set = {t for t in tips}

    tips_under: dict[Node, list[Node]] = {}

    def collect(node: Node) -> list[Node]:
        if node.is_leaf:
            got = [node] if node in tip_set else []
            tips_under[node] = got
            return got
        got: list[Node] = []
        for c in node.children:
            got.extend(collect(c))
        tips_under[node] = got
        return got

    collect(tree_root)

    best_meta: Optional[tuple[int, float]] = None
    best_left: Optional[list[Node]] = None
    best_right: Optional[list[Node]] = None

    def consider(right_tips: list[Node]) -> None:
        nonlocal best_meta, best_left, best_right
        rn = len(right_tips)
        if rn == 0 or rn == n:
            return
        right_set = set(right_tips)
        left_tips = [t for t in tips if t not in right_set]
        # Both children must meet min cluster size (no undersized / noise-only splits).
        if len(left_tips) < min_size or len(right_tips) < min_size:
            return
        d0 = max_depth_of_tips(left_tips, root_dist)
        d1 = max_depth_of_tips(right_tips, root_dist)
        score = max(d0, d1)
        meta = (0, score)
        if _is_better_bip(meta, best_meta):
            best_meta = meta
            best_left = left_tips
            best_right = right_tips

    for _node, under in tips_under.items():
        consider(under)

    # Fallback: partition children of cluster LCA
    if best_meta is None and len(tips) >= 2 * min_size:
        anc = tips[0]
        for t in tips[1:]:
            anc = lca(anc, t)
        kids = list(anc.children) if anc.children else []
        if 2 <= len(kids) <= 20:
            child_lists = [tips_under.get(ch, []) for ch in kids]
            m = len(child_lists)
            limit = (1 << m) - 1
            for mask in range(1, limit):
                right: list[Node] = []
                for b in range(m):
                    if mask & (1 << b):
                        right.extend(child_lists[b])
                consider(right)

    if best_meta is None or best_left is None or best_right is None:
        return None
    return best_left, best_right, best_meta[1]


def find_leaf_to_split(
    root: KDpisLeaf,
    min_size: int,
    require_dpi_above: Optional[int],
) -> Optional[KDpisLeaf]:
    # Need room for two children each ≥ min_size.
    min_split_n = max(2, 2 * min_size)
    leaves = collect_kdpis_leaves(root)
    best: Optional[KDpisLeaf] = None
    best_key = -math.inf
    for leaf in leaves:
        if leaf.unsplittable or leaf.n < min_split_n:
            continue
        if require_dpi_above is not None:
            if leaf.dpi is None or not (leaf.dpi > require_dpi_above):
                continue
        elif leaf.dpi is None:
            continue
        key = float(leaf.dpi) if leaf.dpi is not None else -math.inf
        if key > best_key or (key == best_key and best is not None and leaf.n > best.n):
            best_key = key
            best = leaf
    return best


def split_leaf(
    leaf: KDpisLeaf,
    tree_root: Node,
    root_dist: dict[Node, float],
    dsr: float,
    min_size: int,
) -> bool:
    if leaf.n < 2 * min_size or not leaf.is_atomic or leaf.unsplittable:
        return False
    bip = best_bipartition(leaf.tips, tree_root, root_dist, min_size)
    if bip is None:
        leaf.unsplittable = True
        return False
    left_tips, right_tips, _score = bip
    leaf.unsplittable = False
    leaf.left = make_kdpis_leaf(left_tips, root_dist, dsr)
    leaf.right = make_kdpis_leaf(right_tips, root_dist, dsr)
    return True


def long_branch_noise_tips(tips: list[Node], max_branch: float) -> tuple[list[Node], list[Node]]:
    keep, noise = [], []
    for t in tips:
        if (t.length or 0.0) > max_branch:
            noise.append(t)
        else:
            keep.append(t)
    return keep, noise


def assign_cluster_ids(
    root: KDpisLeaf,
    tip_order: dict[str, int],
    min_size: int,
    pre_noise: list[Node],
) -> dict[str, int]:
    """name → cluster id (−1 = noise). Non-noise numbered 1..k by phylogeny tip order.

    After min-size filtering, every remaining singleton is reclassified as noise
    before tree-order ID assignment.
    """
    leaves = collect_kdpis_leaves(root)
    noise_leaves = [L for L in leaves if L.n < min_size]
    cluster_leaves = [L for L in leaves if L.n >= min_size]

    # Reclassify singletons as noise before tree-order renumbering.
    kept: list[KDpisLeaf] = []
    for L in cluster_leaves:
        if L.n == 1:
            noise_leaves.append(L)
        else:
            kept.append(L)
    cluster_leaves = kept

    def leaf_y(L: KDpisLeaf) -> float:
        ys = [tip_order.get(t.name or "", 10**9) for t in L.tips]
        return min(ys) if ys else 10**9

    cluster_leaves.sort(key=leaf_y)
    name_to_cl: dict[str, int] = {}
    for L in noise_leaves:
        L.cluster_id = -1
        L.is_noise = True
        for t in L.tips:
            if t.name:
                name_to_cl[normalize_tip_name(t.name)] = -1
    for t in pre_noise:
        if t.name:
            name_to_cl[normalize_tip_name(t.name)] = -1
    next_id = 1
    for L in cluster_leaves:
        L.cluster_id = next_id
        L.is_noise = False
        for t in L.tips:
            if t.name:
                name_to_cl[normalize_tip_name(t.name)] = next_id
        next_id += 1
    return name_to_cl


def run_kdpis_auto(
    tree: Tree,
    target_dpi: int,
    dsr: float = DEFAULT_DSR,
    min_size: int = DEFAULT_MIN_SIZE,
    max_splits: int = DEFAULT_MAX_SPLITS,
    remove_long_branches: bool = True,
) -> tuple[dict[str, int], dict[int, Optional[int]], dict[int, int], KDpisLeaf]:
    """
    Returns:
      name→clusterId, clusterId→dpi, clusterId→count, k-DPIs root
    """
    set_parents(tree.root)
    root_dist = node_to_root_distance(tree.root)
    all_tips = [t for t in tree.tips() if t.name]
    etd = 2.0 * dsr * target_dpi
    if remove_long_branches:
        cluster_tips, pre_noise = long_branch_noise_tips(all_tips, etd)
    else:
        cluster_tips, pre_noise = list(all_tips), []
    if len(cluster_tips) < 1:
        cluster_tips, pre_noise = list(all_tips), []

    k_root = make_kdpis_leaf(cluster_tips, root_dist, dsr)
    steps = 0
    while steps < max_splits:
        leaf = find_leaf_to_split(k_root, min_size, require_dpi_above=target_dpi)
        if leaf is None:
            break
        if not split_leaf(leaf, tree.root, root_dist, dsr, min_size):
            continue
        steps += 1

    tip_order = {t.name: i for i, t in enumerate(tree.tips()) if t.name}
    name_to_cl = assign_cluster_ids(k_root, tip_order, min_size, pre_noise)

    # DPI / counts per cluster
    by_cl: dict[int, list[Node]] = defaultdict(list)
    tip_by_name = {normalize_tip_name(t.name): t for t in all_tips if t.name}
    for nm, cid in name_to_cl.items():
        if cid == -1:
            continue
        t = tip_by_name.get(nm)
        if t:
            by_cl[cid].append(t)
    cl_dpi = {cid: dpi_days_from_mean(mean_pairwise_path(ts, root_dist), dsr) for cid, ts in by_cl.items()}
    cl_count = {cid: len(ts) for cid, ts in by_cl.items()}
    return name_to_cl, cl_dpi, cl_count, k_root


# ---------------------------------------------------------------------------
# SVG
# ---------------------------------------------------------------------------

def escape_xml(s: str) -> str:
    return (
        str(s)
        .replace("&", "&amp;")
        .replace("<", "&lt;")
        .replace(">", "&gt;")
        .replace('"', "&quot;")
    )


def cluster_color(cid: int) -> str:
    return COLOR_PALETTE[cid % len(COLOR_PALETTE)]


def group_color(gid: int) -> str:
    return COLOR_PALETTE[(gid + 1) % len(COLOR_PALETTE)]


def format_branch(x: float) -> str:
    if x == 0:
        return "0"
    ax = abs(x)
    if ax >= 1:
        return f"{x:.4g}"
    if ax >= 0.001:
        return f"{x:.5f}".rstrip("0").rstrip(".")
    return f"{x:.3e}"


def pick_phylo_scale_bar_value(max_depth: float) -> Optional[float]:
    """Nice scale-bar length ≈ max_depth/3 (aliViz pickPhyloScaleBarValue)."""
    if max_depth is None or not math.isfinite(max_depth) or max_depth <= 0:
        return None
    target = max_depth / 3.0
    log10 = math.log10(target)
    if not math.isfinite(log10):
        return None
    exp = math.floor(log10)
    base = 10.0 ** exp
    f = target / base
    if f <= 1:
        nf = 1
    elif f <= 2:
        nf = 2
    elif f <= 5:
        nf = 5
    else:
        nf = 10
    v = nf * base
    if v > max_depth:
        v = base
        while v > max_depth and v > 0:
            v /= 10.0
        if v <= 0:
            v = max_depth
    return v


def compute_svg_x_scale_dpi(dpi_days: int, dsr: float) -> tuple[float, float, str]:
    """
    Legacy DPI scale (aliViz --scale dpi): map expected = DPI×dsr to DPI_MARK_PX (520/16).

    Returns (x_scale, scale_bar_val, scale_bar_label).
    Kept so callers can revert with --scale dpi.
    """
    expected = dpi_days * dsr
    if expected > 0:
        x_scale = DPI_MARK_PX / expected
    else:
        x_scale = DPI_MARK_PX / (DEFAULT_DPI_DAYS * DEFAULT_DSR)
        expected = DEFAULT_DPI_DAYS * DEFAULT_DSR
    label = f"({dpi_days} DPI @ {dsr} = {format_branch(expected)})"
    return x_scale, expected, label


def compute_svg_x_scale_auto(
    max_depth: float,
    dpi_days: Optional[int],
    dsr: float,
) -> tuple[float, Optional[float], Optional[str], bool]:
    """
    AUTO scale (aliViz default): fit tree into BRANCH_W; if dpi_days is set, also fit
    expected DPI depth so the vertical mark stays visible, and use DPI scale-bar label.

    Returns (x_scale, scale_bar_val, scale_bar_label, draw_dpi_mark).
    """
    max_d = max_depth if max_depth > 0 else 1.0
    expected = (dpi_days * dsr) if (dpi_days is not None and dpi_days > 0 and dsr > 0) else None
    if expected is not None and expected > 0:
        fit_depth = max(max_d, expected)
        x_scale = BRANCH_W / fit_depth
        label = f"({dpi_days} DPI @ {dsr} = {format_branch(expected)})"
        return x_scale, expected, label, True
    x_scale = BRANCH_W / max_d
    bar = pick_phylo_scale_bar_value(max_d)
    return x_scale, bar, None, False


def build_svg(
    tree: Tree,
    max_depth: float,
    dpi_days: int,
    dsr: float,
    title: str,
    name_to_gid: dict[str, int],
    gid_to_label: dict[int, str],
    group_counts: dict[int, int],
    group_dpi: dict[int, Optional[int]],
    name_to_cl: dict[str, int],
    cl_dpi: dict[int, Optional[int]],
    cl_count: dict[int, int],
    founder_name: Optional[str],
    scale_mode: str = DEFAULT_SCALE_MODE,
    dpi_detected: bool = True,
) -> str:
    """
    Linear tree SVG.

    scale_mode:
      - "auto" (default): fit tree (+ DPI mark when dpi_detected) into BRANCH_W
      - "dpi": legacy fixed scale with expected depth at BRANCH_W/16
    """
    mode = (scale_mode or DEFAULT_SCALE_MODE).lower()
    if mode not in ("auto", "dpi"):
        mode = DEFAULT_SCALE_MODE

    # --- x-scale / scale bar ---
    # Legacy DPI path retained for --scale dpi (maps DPI×dsr → 520/16 px).
    if mode == "dpi":
        x_scale, scale_bar_val, scale_bar_label = compute_svg_x_scale_dpi(dpi_days, dsr)
        draw_dpi_mark = scale_bar_val is not None and scale_bar_val > 0
    else:
        overlay_days = dpi_days if dpi_detected else None
        x_scale, scale_bar_val, scale_bar_label, draw_dpi_mark = compute_svg_x_scale_auto(
            max_depth, overlay_days, dsr
        )

    padding = 15.0
    left_pad = 5.0
    tips = tree.tips()
    n_rows = max(1, len(tips))
    # Extra top row for the scale bar (aliViz parks it on the REF row above the first tip).
    scale_bar_rows = 1
    svg_h = (n_rows + scale_bar_rows) * ROW_H
    label_x0 = left_pad + 5 + BRANCH_W + 14
    max_label_w = 420.0
    title_space = 44.0
    plot_w = label_x0 + max_label_w + padding
    total_w = plot_w * (8.0 / 7.0)
    legend_w = total_w - plot_w
    total_h = svg_h + padding * 2 + title_space
    clip_right = left_pad + BRANCH_W

    def nx(n: Node) -> float:
        return left_pad + n.x_depth * x_scale

    def ny(n: Node) -> float:
        # Shift tips down one row so row 0 is free for the scale bar.
        return (n.y_row + scale_bar_rows) * ROW_H + ROW_H / 2.0

    def tip_fill(name: Optional[str]) -> str:
        if not name:
            return "#374151"
        cid = name_to_cl.get(normalize_tip_name(name))
        if cid is None:
            return "#374151"
        if cid == -1:
            return NOISE_COLOR
        return cluster_color(cid)

    def label_fill(name: Optional[str]) -> str:
        if not name:
            return "#374151"
        nm = normalize_tip_name(name)
        if founder_name and nm == normalize_tip_name(founder_name):
            return "#ff00ff"
        gid = name_to_gid.get(nm)
        if gid is None:
            return "#374151"
        return group_color(gid)

    def display_label(name: Optional[str]) -> str:
        if not name:
            return ""
        nm = normalize_tip_name(name)
        lab = name
        if founder_name and nm == normalize_tip_name(founder_name):
            lab = f"{name} [Founder]"
        return escape_xml(lab)

    num_groups = len(gid_to_label)
    num_clusters = sum(1 for c in set(name_to_cl.values()) if c != -1)
    noise_n = sum(1 for c in name_to_cl.values() if c == -1)
    max_d_str = format_branch(max_depth)
    subtitle = (
        f"Infer: loaded tree | s={len(tips)} | Max depth: {max_d_str} | "
        f"Cluster: k-DPIs | founder={escape_xml(founder_name or '—')}"
    )

    parts: list[str] = []
    parts.append('<?xml version="1.0" encoding="UTF-8"?>\n')
    parts.append(
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{total_w}" height="{total_h}" '
        f'viewBox="0 0 {total_w} {total_h}">\n'
    )
    parts.append(f'<rect width="{total_w}" height="{total_h}" fill="white"/>\n')
    parts.append('<g font-family="sans-serif" fill="#111827" text-anchor="start">\n')
    parts.append(f'<text x="{padding + 5}" y="20" font-size="13" font-weight="700">{escape_xml(title)}</text>\n')
    parts.append(f'<text x="{padding + 5}" y="36" font-size="11">{subtitle}</text>\n')
    parts.append("</g>\n")

    parts.append(
        f'<g transform="translate({padding},{padding + title_space})" '
        f'stroke="#374151" stroke-width="1" fill="none">\n'
    )

    def draw_horiz(x1: float, y: float, x2: float) -> None:
        if x2 <= clip_right:
            parts.append(f'<line x1="{x1}" y1="{y}" x2="{x2}" y2="{y}"/>\n')
        elif x1 < clip_right:
            parts.append(f'<line x1="{x1}" y1="{y}" x2="{clip_right}" y2="{y}"/>\n')

    internals: list[Node] = []
    leaves_draw: list[Node] = []

    def draw_node(n: Node) -> None:
        x, y = nx(n), ny(n)
        if not n.is_leaf:
            min_y, max_y = math.inf, -math.inf
            for c in n.children:
                draw_node(c)
                cx, cy = nx(c), ny(c)
                draw_horiz(x, cy, cx)
                min_y = min(min_y, cy)
                max_y = max(max_y, cy)
            if min_y != math.inf and x <= clip_right:
                parts.append(f'<line x1="{x}" y1="{min_y}" x2="{x}" y2="{max_y}"/>\n')
            internals.append(n)
        else:
            leaves_draw.append(n)

    draw_node(tree.root)

    if draw_dpi_mark and scale_bar_val is not None and scale_bar_val > 0:
        expected_x = left_pad + scale_bar_val * x_scale
        parts.append(
            f'<line x1="{expected_x}" y1="0" x2="{expected_x}" y2="{svg_h}" '
            f'stroke="#6b7280" stroke-width="1" stroke-dasharray="3 3"/>\n'
        )

    parts.append('<g stroke="#374151" fill="none">\n')
    for n in internals:
        x, y = nx(n), ny(n)
        if x > clip_right:
            continue
        parts.append(f'<circle cx="{x}" cy="{y}" r="3" fill="#374151" stroke="none"/>\n')
    d = 10.0
    for n in leaves_draw:
        x, y = nx(n), ny(n)
        fill = tip_fill(n.name)
        if x > clip_right:
            cx = clip_right - d
            parts.append(
                f'<polygon points="{cx},{y - d} {clip_right},{y} {cx},{y + d}" '
                f'fill="{fill}" stroke="#374151" stroke-width="0.5"/>\n'
            )
            continue
        pts = f"{x},{y - d} {x + d},{y} {x},{y + d} {x - d},{y}"
        parts.append(f'<polygon points="{pts}" fill="{fill}" stroke="#374151" stroke-width="0.5"/>\n')
    parts.append("</g>\n")

    parts.append('<g fill="none" stroke="#374151" stroke-width="1">\n')
    for n in leaves_draw:
        x, y = nx(n), ny(n)
        x_start = min(x, left_pad + BRANCH_W)
        stroke = label_fill(n.name)
        parts.append(
            f'<line x1="{x_start}" y1="{y}" x2="{label_x0 - 6}" y2="{y}" '
            f'stroke-dasharray="2 4" stroke="{stroke}"/>\n'
        )
    parts.append("</g>\n")
    parts.append("</g>\n")

    # Labels
    parts.append(
        f'<g transform="translate({padding},{padding + title_space})" '
        f'font-family="sans-serif" font-size="11">\n'
    )

    def draw_labels(n: Node) -> None:
        if n.is_leaf:
            if n.name:
                y = ny(n)
                parts.append(
                    f'<text x="{label_x0}" y="{y}" fill="{label_fill(n.name)}" '
                    f'font-weight="600" dominant-baseline="middle">{display_label(n.name)}</text>\n'
                )
            return
        for c in n.children:
            draw_labels(c)

    draw_labels(tree.root)
    parts.append("</g>\n")

    # Scale bar (centred in the top spacer row, above the first tip)
    if scale_bar_val is not None and scale_bar_val > 0:
        parts.append(f'<g transform="translate({padding},{padding + title_space})">\n')
        px_len = scale_bar_val * x_scale
        x0, y0 = 5.0, ROW_H / 2.0
        label = scale_bar_label if scale_bar_label else format_branch(scale_bar_val)
        parts.append('<g fill="none" stroke="#111827" stroke-width="2">\n')
        parts.append(f'<line x1="{x0}" y1="{y0}" x2="{x0 + px_len}" y2="{y0}"/>\n')
        parts.append(f'<line x1="{x0}" y1="{y0 - 4}" x2="{x0}" y2="{y0 + 4}"/>\n')
        parts.append(f'<line x1="{x0 + px_len}" y1="{y0 - 4}" x2="{x0 + px_len}" y2="{y0 + 4}"/>\n')
        parts.append("</g>\n")
        parts.append(
            f'<text x="{x0 + px_len + 8}" y="{y0}" font-family="sans-serif" '
            f'font-size="11" fill="#111827" dominant-baseline="middle">{escape_xml(label)}</text>\n'
        )
        parts.append("</g>\n")

    # Legend
    legend_x = plot_w
    legend_y = title_space
    pad = 8.0
    sw = 14.0
    row_h = 20.0
    y = legend_y + pad + 14
    x_leg = legend_x + pad
    leg: list[str] = []
    leg.append(
        f'<text x="{x_leg}" y="{y}" font-family="sans-serif" font-size="14" '
        f'font-weight="700" fill="#374151">Legend</text>\n'
    )
    y += 22
    leg.append('<g font-family="sans-serif" font-size="11">\n')
    if gid_to_label:
        leg.append(f'<text x="{x_leg}" y="{y}" font-weight="600" fill="#374151">Groups:</text>\n')
        y += 16
        for gid in sorted(gid_to_label.keys(), key=lambda g: gid_to_label[g]):
            color = group_color(gid)
            n = group_counts.get(gid, 0)
            dpi = group_dpi.get(gid)
            lab = f"{gid_to_label[gid]} ({n}, {dpi})" if dpi is not None else f"{gid_to_label[gid]} ({n})"
            leg.append(
                f'<rect x="{x_leg}" y="{y - sw / 2 + 1}" width="{sw}" height="{sw}" '
                f'fill="{color}" stroke="#d1d5db" stroke-width="1" rx="2"/>\n'
            )
            leg.append(f'<text x="{x_leg + sw + 8}" y="{y + 4}" fill="#374151">{escape_xml(lab)}</text>\n')
            y += row_h
        y += 6
    cluster_ids = sorted(c for c in set(name_to_cl.values()) if c != -1)
    if cluster_ids or noise_n:
        leg.append(f'<text x="{x_leg}" y="{y}" font-weight="600" fill="#374151">Clusters:</text>\n')
        y += 16
        for cid in cluster_ids:
            color = cluster_color(cid)
            n = cl_count.get(cid, 0)
            dpi = cl_dpi.get(cid)
            lab = f"{cid} ({n}, {dpi})" if dpi is not None else f"{cid} ({n})"
            leg.append(
                f'<rect x="{x_leg}" y="{y - sw / 2 + 1}" width="{sw}" height="{sw}" '
                f'fill="{color}" stroke="#d1d5db" stroke-width="1" rx="2"/>\n'
            )
            leg.append(f'<text x="{x_leg + sw + 8}" y="{y + 4}" fill="#374151">{escape_xml(lab)}</text>\n')
            y += row_h
        if noise_n:
            leg.append(
                f'<rect x="{x_leg}" y="{y - sw / 2 + 1}" width="{sw}" height="{sw}" '
                f'fill="{NOISE_COLOR}" stroke="#d1d5db" stroke-width="1" rx="2"/>\n'
            )
            leg.append(
                f'<text x="{x_leg + sw + 8}" y="{y + 4}" fill="#374151">'
                f'{escape_xml(f"Noise ({noise_n})")}</text>\n'
            )
            y += row_h
    leg.append("</g>\n")
    legend_h = max(y - legend_y + pad, 24)
    parts.append(
        f'<rect x="{legend_x}" y="{legend_y}" width="{legend_w}" height="{legend_h}" '
        f'fill="#fafafa" stroke="#e5e7eb" stroke-width="1"/>\n'
    )
    parts.extend(leg)
    parts.append("</svg>")
    _ = (num_groups, num_clusters)
    return "".join(parts)


# ---------------------------------------------------------------------------
# Pipeline per tree
# ---------------------------------------------------------------------------

def process_tree_file(
    path: Path,
    out_dir: Path,
    delimiter: str,
    field: int,
    dsr: float,
    min_size: int,
    max_splits: int,
    remove_long_branches: bool,
    scale_mode: str = DEFAULT_SCALE_MODE,
) -> tuple[Path, dict[str, object]]:
    text = path.read_text(encoding="utf-8", errors="replace")
    tree = parse_newick(text)
    set_parents(tree.root)

    dpi_detected = detect_dpi_from_filename(path.name) is not None
    dpi_days = parse_dpi_from_filename(path.name)
    tip_names = [t.name for t in tree.tips() if t.name]
    name_to_gid, gid_to_label, members = group_tips(tip_names, delimiter, field)
    if not members:
        raise ValueError(f"No tips could be grouped with delimiter={delimiter!r} field={field}")

    # Founder = tree-path medoid of the smallest group
    sgid = smallest_group_id(members)
    assert sgid is not None
    tip_by_name = {normalize_tip_name(t.name): t for t in tree.tips() if t.name}
    group_nodes = [tip_by_name[nm] for nm in members[sgid] if nm in tip_by_name]
    if not group_nodes:
        raise ValueError("Smallest group has no matching tree tips")
    root_dist = node_to_root_distance(tree.root)
    founder = tree_medoid(group_nodes, root_dist)
    founder_name = founder.name

    tree = reroot_on_leaf(tree, founder)
    # Re-find founder tip after reroot (object identity preserved if same Node)
    ladderize_by_depth(tree)
    max_depth = layout_tree(tree)
    set_parents(tree.root)
    root_dist = node_to_root_distance(tree.root)

    # Group DPI after reroot/ladderize
    tip_by_name = {normalize_tip_name(t.name): t for t in tree.tips() if t.name}
    group_counts = {gid: len(mems) for gid, mems in members.items()}
    group_dpi: dict[int, Optional[int]] = {}
    for gid, mems in members.items():
        nodes = [tip_by_name[nm] for nm in mems if nm in tip_by_name]
        group_dpi[gid] = dpi_days_from_mean(mean_pairwise_path(nodes, root_dist), dsr)

    name_to_cl, cl_dpi, cl_count, _kroot = run_kdpis_auto(
        tree,
        target_dpi=dpi_days,
        dsr=dsr,
        min_size=min_size,
        max_splits=max_splits,
        remove_long_branches=remove_long_branches,
    )
    # Re-layout after clustering (topology unchanged; y_row already set)
    max_depth = layout_tree(tree)

    num_groups = len(gid_to_label)
    num_clusters = len([c for c in set(name_to_cl.values()) if c != -1])
    base = path.stem
    out_name = f"{base}_gr-{num_groups}_cl-{num_clusters}.svg"
    out_path = out_dir / out_name

    svg = build_svg(
        tree=tree,
        max_depth=max_depth,
        dpi_days=dpi_days,
        dsr=dsr,
        title=path.name,
        name_to_gid=name_to_gid,
        gid_to_label=gid_to_label,
        group_counts=group_counts,
        group_dpi=group_dpi,
        name_to_cl=name_to_cl,
        cl_dpi=cl_dpi,
        cl_count=cl_count,
        founder_name=founder_name,
        scale_mode=scale_mode,
        dpi_detected=dpi_detected,
    )
    out_path.write_text(svg, encoding="utf-8")
    record = {
        "CAPid": parse_capid_from_filename(path.name),
        "Groups": num_groups,
        "Clusters": num_clusters,
    }
    return out_path, record


def write_housekeeping_csv(out_dir: Path, records: list[dict[str, object]]) -> Path:
    """Write housekeeping.csv sorted by CAPid (numeric)."""
    rows = sorted(records, key=lambda r: capid_sort_key(str(r.get("CAPid") or "")))
    path = out_dir / "housekeeping.csv"
    with path.open("w", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=["CAPid", "Groups", "Clusters"])
        writer.writeheader()
        for row in rows:
            writer.writerow({
                "CAPid": row.get("CAPid") or "",
                "Groups": row.get("Groups", 0),
                "Clusters": row.get("Clusters", 0),
            })
    return path


def _pie_slice_path(cx: float, cy: float, r: float, a0: float, a1: float) -> str:
    """SVG path for a pie slice from angle a0 to a1 (radians, 0 = east, CCW)."""
    if abs(a1 - a0) >= 2 * math.pi - 1e-9:
        # Full circle as two semicircles
        x1 = cx + r
        y1 = cy
        return (
            f"M {cx} {cy} L {x1} {y1} "
            f"A {r} {r} 0 1 0 {cx - r} {cy} "
            f"A {r} {r} 0 1 0 {x1} {y1} Z"
        )
    x0 = cx + r * math.cos(a0)
    y0 = cy + r * math.sin(a0)
    x1 = cx + r * math.cos(a1)
    y1 = cy + r * math.sin(a1)
    large = 1 if (a1 - a0) > math.pi else 0
    return f"M {cx} {cy} L {x0} {y0} A {r} {r} 0 {large} 1 {x1} {y1} Z"


def cluster_count_bin_label(n_clusters: int) -> str:
    """
    Bin label for cluster count: 1, 2, 3–4, 5–8, 9–16, … (powers of two).
    Counts ≤0 are binned with 1.
    """
    n = max(1, int(n_clusters))
    if n == 1:
        return "1"
    if n == 2:
        return "2"
    hi = 4
    while n > hi:
        hi *= 2
    lo = hi // 2 + 1
    return f"{lo}-{hi}"


def cluster_count_bin_sort_key(label: str) -> int:
    """Order bins: 1, 2, 3-4, 5-8, …"""
    m = re.match(r"^(\d+)-", label)
    if m:
        return int(m.group(1))
    m2 = re.match(r"^(\d+)$", label)
    return int(m2.group(1)) if m2 else 10**9


def write_housekeeping_pie_svg(out_dir: Path, records: list[dict[str, object]]) -> Path:
    """
    Write housekeeping.svg: pie chart of CAP samples by cluster-count bins
    (1, 2, 3–4, 5–8, 9–16, …).
    """
    counts: dict[str, int] = defaultdict(int)
    for row in records:
        try:
            n_cl = int(row.get("Clusters") or 0)
        except (TypeError, ValueError):
            n_cl = 0
        counts[cluster_count_bin_label(n_cl)] += 1

    # Only non-empty bins, ordered
    bins = sorted(counts.keys(), key=cluster_count_bin_sort_key)
    slices = [(lab, counts[lab]) for lab in bins if counts[lab] > 0]
    total = sum(c for _, c in slices)

    palette = [
        "#3b82f6",
        "#10b981",
        "#f59e0b",
        "#8b5cf6",
        "#06b6d4",
        "#f97316",
        "#6366f1",
        "#14b8a6",
        "#a855f7",
        "#22c55e",
        "#64748b",
    ]

    def bin_color(idx: int) -> str:
        return palette[idx % len(palette)]

    legend_rows = max(len(slices), 1)
    w = 520.0
    h = max(360.0, 80.0 + legend_rows * 24.0 + 40.0)
    cx, cy, r = 180.0, max(180.0, h / 2.0), 120.0

    parts: list[str] = []
    parts.append('<?xml version="1.0" encoding="UTF-8"?>\n')
    parts.append(
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{w}" height="{h}" '
        f'viewBox="0 0 {w} {h}">\n'
    )
    parts.append(f'<rect width="{w}" height="{h}" fill="white"/>\n')
    parts.append(
        '<text x="24" y="28" font-family="sans-serif" font-size="16" font-weight="700" '
        'fill="#111827">CAP samples by cluster-count bin</text>\n'
    )

    if total == 0:
        parts.append(
            '<text x="24" y="80" font-family="sans-serif" font-size="13" fill="#6b7280">'
            "No CAP records to plot.</text>\n"
        )
    elif len(slices) == 1:
        color = bin_color(0)
        parts.append(
            f'<circle cx="{cx}" cy="{cy}" r="{r}" fill="{color}" '
            f'stroke="#ffffff" stroke-width="2"/>\n'
        )
        parts.append(
            f'<text x="{cx}" y="{cy}" text-anchor="middle" dominant-baseline="middle" '
            f'font-family="sans-serif" font-size="18" font-weight="700" fill="#ffffff">'
            f"100%</text>\n"
        )
    else:
        start = -math.pi / 2
        angle = start
        for i, (lab, cnt) in enumerate(slices):
            frac = cnt / total
            a0 = angle
            a1 = angle + frac * 2 * math.pi
            color = bin_color(i)
            parts.append(
                f'<path d="{_pie_slice_path(cx, cy, r, a0, a1)}" '
                f'fill="{color}" stroke="#ffffff" stroke-width="2"/>\n'
            )
            if frac >= 0.04:
                mid = (a0 + a1) / 2.0
                lx = cx + (r * 0.55) * math.cos(mid)
                ly = cy + (r * 0.55) * math.sin(mid)
                pct = 100.0 * frac
                parts.append(
                    f'<text x="{lx:.1f}" y="{ly:.1f}" text-anchor="middle" '
                    f'dominant-baseline="middle" font-family="sans-serif" font-size="12" '
                    f'font-weight="700" fill="#ffffff">{pct:.0f}%</text>\n'
                )
            angle = a1

    # Legend
    lx0, ly0 = 340.0, 70.0
    sw = 16.0
    parts.append('<g font-family="sans-serif" font-size="13" fill="#374151">\n')
    y = ly0
    for i, (lab, cnt) in enumerate(slices):
        color = bin_color(i)
        pct = (100.0 * cnt / total) if total else 0.0
        if "-" in lab:
            desc = f"Clusters {lab}: {cnt} ({pct:.0f}%)"
        else:
            desc = f"Clusters = {lab}: {cnt} ({pct:.0f}%)"
        parts.append(
            f'<rect x="{lx0}" y="{y}" width="{sw}" height="{sw}" fill="{color}" '
            f'stroke="#d1d5db" rx="2"/>\n'
        )
        parts.append(f'<text x="{lx0 + sw + 8}" y="{y + 13}">{escape_xml(desc)}</text>\n')
        y += 24.0
    parts.append(
        f'<text x="{lx0}" y="{y + 16}" font-size="12" fill="#6b7280">'
        f"n = {total} CAP sample(s)</text>\n"
    )
    parts.append("</g>\n")
    parts.append("</svg>\n")

    path = out_dir / "housekeeping.svg"
    path.write_text("".join(parts), encoding="utf-8")
    return path


def iter_tree_files(input_dir: Path) -> list[Path]:
    files = []
    for p in sorted(input_dir.iterdir()):
        if p.is_file() and p.suffix.lower() in TREE_EXTENSIONS:
            files.append(p)
    return files


def main(argv: Optional[list[str]] = None) -> int:
    ap = argparse.ArgumentParser(
        description="Batch k-DPIs clustering + SVG export (aliViz-compatible)."
    )
    ap.add_argument("input_dir", type=Path, help="Directory of Newick tree files")
    ap.add_argument("output_dir", type=Path, help="Directory for output SVG files")
    ap.add_argument("--delimiter", default=DEFAULT_DELIMITER, help="Group name delimiter (default _)")
    ap.add_argument(
        "--field",
        type=int,
        default=DEFAULT_FIELD,
        help="1-based name field for grouping (default 3, same as aliViz)",
    )
    ap.add_argument("--dsr", type=float, default=DEFAULT_DSR, help=f"Daily substitution rate (default {DEFAULT_DSR})")
    ap.add_argument("--min-size", type=int, default=DEFAULT_MIN_SIZE, help=f"Min cluster size (default {DEFAULT_MIN_SIZE})")
    ap.add_argument("--max-splits", type=int, default=DEFAULT_MAX_SPLITS, help=f"Max k-DPIs splits (default {DEFAULT_MAX_SPLITS})")
    ap.add_argument(
        "--scale",
        choices=("auto", "dpi"),
        default=DEFAULT_SCALE_MODE,
        help="SVG branch scale: auto (default, fit tree + DPI mark) or dpi (legacy 1/16 mark)",
    )
    ap.add_argument(
        "--no-long-branch-noise",
        action="store_true",
        help="Disable removal of tips with incoming branch > etd",
    )
    args = ap.parse_args(argv)

    in_dir: Path = args.input_dir
    out_dir: Path = args.output_dir
    if not in_dir.is_dir():
        print(f"Error: input_dir is not a directory: {in_dir}", file=sys.stderr)
        return 1
    out_dir.mkdir(parents=True, exist_ok=True)

    files = iter_tree_files(in_dir)
    if not files:
        print(f"No tree files found in {in_dir} (extensions: {sorted(TREE_EXTENSIONS)})", file=sys.stderr)
        return 1

    ok, fail = 0, 0
    records: list[dict[str, object]] = []
    for path in files:
        try:
            out, record = process_tree_file(
                path,
                out_dir,
                delimiter=args.delimiter,
                field=args.field,
                dsr=args.dsr,
                min_size=args.min_size,
                max_splits=args.max_splits,
                remove_long_branches=not args.no_long_branch_noise,
                scale_mode=args.scale,
            )
            records.append(record)
            print(f"OK  {path.name} -> {out.name}  ({record['CAPid'] or '—'} gr={record['Groups']} cl={record['Clusters']})")
            ok += 1
        except Exception as exc:
            print(f"FAIL {path.name}: {exc}", file=sys.stderr)
            fail += 1

    hk = write_housekeeping_csv(out_dir, records)
    pie = write_housekeeping_pie_svg(out_dir, records)
    print(f"Housekeeping: {hk.name} ({len(records)} record(s)), {pie.name}")
    print(f"Done: {ok} written, {fail} failed → {out_dir}")
    return 0 if fail == 0 else 2


if __name__ == "__main__":
    raise SystemExit(main())
