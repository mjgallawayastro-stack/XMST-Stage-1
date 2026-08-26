#!/usr/bin/env python3
"""
Nested XMST distance-uncertainty experiment for Paper I.

Reproduces the authoritative empirical-template Monte Carlo geometry and then,
for each nominal field, perturbs stellar distances repeatedly using Q25's
published asymmetric 16th/50th/84th-percentile distance intervals.

Important design choice for moved template stars:
    the source Q25 star's *fractional* lower/upper distance uncertainty is
    preserved after the template is rotated/translated to its nominal MC
    position. Background stars retain their original Q25 fractional uncertainty.

Every uncertainty realisation independently recomputes:
    Cartesian positions -> Delaunay graph -> Euclidean MST -> percolation limit
    -> two-class Jenks break -> fracture -> retained groups (N >= min group size)

The uncertainty experiment is spatial only: no A_V branch pruning is applied.
"""
from __future__ import annotations

import argparse
import json
import math
import os
import platform
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd
import scipy
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components, minimum_spanning_tree
from scipy.spatial import Delaunay
from scipy.spatial.transform import Rotation
from scipy.stats import norm

BASE = Path(__file__).resolve().parent
DEFAULT_INPUT = BASE / "input"
DEFAULT_OUTPUT = BASE / "output"

OUTER_SEED = 20260725
INNER_SEED = 20260818

CFG = {
    "min_group_size": 10,
    "overlap_distance_min_pc": 57.0,
    "overlap_distance_max_pc": 83.0,
    "centre_min_separation_pc": 120.0,
    "centre_distance_min_pc": 250.0,
    "centre_distance_max_pc": 900.0,
}


@dataclass
class LoadedData:
    full_df: pd.DataFrame
    field_df: pd.DataFrame
    templates: list[dict]
    template_labels: list[str]


def exact_mst(xyz: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray, int]:
    """Exact 3-D Euclidean MST via the Delaunay graph, matching the MC1000 code."""
    tri = Delaunay(xyz, qhull_options="Qbb Qc Qz Q12")
    s = tri.simplices
    pairs = np.vstack([
        s[:, [0, 1]], s[:, [0, 2]], s[:, [0, 3]],
        s[:, [1, 2]], s[:, [1, 3]], s[:, [2, 3]],
    ])
    pairs = np.sort(pairs, axis=1)
    pairs = np.unique(pairs, axis=0)
    n = len(xyz)
    # Qhull's point-at-infinity can appear as index n under Qz.
    pairs = pairs[(pairs[:, 0] < n) & (pairs[:, 1] < n)]
    w = np.linalg.norm(xyz[pairs[:, 0]] - xyz[pairs[:, 1]], axis=1)
    g = coo_matrix((w, (pairs[:, 0], pairs[:, 1])), shape=(n, n))
    g = g + g.T
    mst = minimum_spanning_tree(g).tocoo()
    u = mst.row.astype(np.int32)
    v = mst.col.astype(np.int32)
    ew = mst.data.astype(float)
    if len(ew) != n - 1:
        raise RuntimeError(f"MST incomplete: {len(ew)} != {n - 1}")
    order = np.argsort(ew, kind="stable")
    return u[order], v[order], ew[order], len(pairs)


def percolation_cut(
    u: np.ndarray, v: np.ndarray, w: np.ndarray, n: int
) -> tuple[float, float, int, int, int]:
    """Finite-size percolation transition used by the authoritative MC1000 run."""
    parent = np.arange(n, dtype=np.int32)
    size = np.ones(n, dtype=np.int32)

    def find(a: int) -> int:
        while parent[a] != a:
            parent[a] = parent[parent[a]]
            a = int(parent[a])
        return a

    largest = 1
    best_jump = -1
    best_i = 0
    best_before = 1
    best_after = 1
    for k in range(len(w)):
        a, b = find(int(u[k])), find(int(v[k]))
        if a == b:
            continue
        if size[a] < size[b]:
            a, b = b, a
        old = largest
        parent[b] = a
        size[a] += size[b]
        largest = max(largest, int(size[a]))
        jump = largest - old
        if jump > best_jump:
            best_jump = jump
            best_i = k
            best_before = old
            best_after = largest
    critical = float(w[best_i])
    previous = float(w[max(0, best_i - 1)])
    operational = 0.5 * (critical + previous)
    return operational, critical, int(best_jump), int(best_before), int(best_after)


def jenks_break(values: np.ndarray) -> tuple[float, float, int]:
    """Exact two-class Jenks split, matching the authoritative MC1000 code."""
    x = np.sort(np.asarray(values, float))
    n = len(x)
    if n < 4:
        return float(np.median(x)), float("nan"), max(1, n // 2)
    cs = np.cumsum(x)
    cs2 = np.cumsum(x * x)
    idx = np.arange(1, n)
    ln = idx
    rn = n - idx
    lsse = cs2[idx - 1] - cs[idx - 1] ** 2 / ln
    rs = cs[-1] - cs[idx - 1]
    rs2 = cs2[-1] - cs2[idx - 1]
    rsse = rs2 - rs ** 2 / rn
    valid = (ln >= 2) & (rn >= 2)
    total = np.where(valid, lsse + rsse, np.inf)
    j = int(np.argmin(total)) + 1
    threshold = float((x[j - 1] + x[j]) / 2.0)
    grand_sse = float(np.sum((x - np.mean(x)) ** 2))
    gvf = 1.0 - float(total[j - 1]) / grand_sse if grand_sse > 0 else float("nan")
    return threshold, gvf, j


def components_from_cut(
    n: int,
    u: np.ndarray,
    v: np.ndarray,
    w: np.ndarray,
    threshold: float,
    min_size: int,
) -> tuple[list[np.ndarray], np.ndarray]:
    keep = w <= threshold
    g = coo_matrix((np.ones(int(keep.sum())), (u[keep], v[keep])), shape=(n, n))
    g = g + g.T
    _, labels = connected_components(g.tocsr(), directed=False)
    counts = np.bincount(labels)
    groups = [np.flatnonzero(labels == gid) for gid, c in enumerate(counts) if c >= min_size]
    # Deterministic output ordering: largest first, then first node.
    groups.sort(key=lambda arr: (-len(arr), int(arr[0])))
    retained_id = np.full(n, -1, dtype=np.int32)
    for gid, members in enumerate(groups):
        retained_id[members] = gid
    return groups, retained_id


def run_xmst(xyz: np.ndarray, min_group_size: int) -> dict:
    u, v, w, candidate_edges = exact_mst(xyz)
    pcut, critical, jump, before, after = percolation_cut(u, v, w, len(xyz))
    sub = w[w <= pcut]
    jcut, gvf, split = jenks_break(sub)
    groups, retained_id = components_from_cut(len(xyz), u, v, w, jcut, min_group_size)
    return {
        "u": u,
        "v": v,
        "w": w,
        "candidate_edges": candidate_edges,
        "percolation_cut_pc": pcut,
        "percolation_critical_edge_pc": critical,
        "percolation_jump": jump,
        "largest_before": before,
        "largest_after": after,
        "subcritical_edge_count": len(sub),
        "jenks_scale_pc": jcut,
        "jenks_gvf": gvf,
        "jenks_split_index": split,
        "groups": groups,
        "group_id": retained_id,
    }


def sample_anchor_indices(field_xyz: np.ndarray, field_dist: np.ndarray, rng: np.random.Generator) -> np.ndarray:
    good = np.flatnonzero(
        (field_dist >= CFG["centre_distance_min_pc"])
        & (field_dist <= CFG["centre_distance_max_pc"])
    )
    if len(good) == 0:
        raise RuntimeError("No anchor candidates")
    return good


def place_template(rel: np.ndarray, centre: np.ndarray, rotation: np.ndarray) -> np.ndarray:
    return rel @ rotation.T + centre


def load_data(input_dir: Path) -> LoadedData:
    q25_path = input_dir / "quintana_2025_ob_stars_full.csv"
    membership_path = input_dir / "q25_group_membership.csv"
    selected_path = input_dir / "selected_template_groups.csv"
    for path in (q25_path, membership_path, selected_path):
        if not path.exists():
            raise FileNotFoundError(path)

    df = pd.read_csv(q25_path, low_memory=False).copy()
    required = {"X", "Y", "Z", "AV", "Dist", "b_Dist", "B_Dist"}
    missing = required.difference(df.columns)
    if missing:
        raise ValueError(f"Q25 input lacks required columns: {sorted(missing)}")
    df["q25_row"] = np.arange(len(df), dtype=np.int32)

    d16 = pd.to_numeric(df["b_Dist"], errors="coerce").to_numpy(float)
    d50 = pd.to_numeric(df["Dist"], errors="coerce").to_numpy(float)
    d84 = pd.to_numeric(df["B_Dist"], errors="coerce").to_numpy(float)
    bad = (~np.isfinite(d16)) | (~np.isfinite(d50)) | (~np.isfinite(d84)) | (d16 <= 0) | (d50 <= 0) | (d84 <= 0) | (d16 > d50) | (d50 > d84)
    if bad.any():
        raise ValueError(f"{int(bad.sum())} Q25 rows have invalid distance percentiles")

    mem = pd.read_csv(membership_path)
    selected = pd.read_csv(selected_path)
    labels: list[str] = []
    templates: list[dict] = []
    selected_rows: list[int] = []

    for rank, r in selected.reset_index(drop=True).iterrows():
        gid = int(r.gid)
        idx = mem.loc[mem.group_label == gid, "row"].to_numpy(int)
        tdf = df.iloc[idx].copy()
        xyz = tdf[["X", "Y", "Z"]].to_numpy(float)
        centre = np.median(xyz, axis=0)
        rel = xyz - centre
        label = f"empirical_template_{rank + 1:02d}"
        labels.append(label)
        templates.append({
            "label": label,
            "gid": gid,
            "indices": idx,
            "q25_rows": tdf["q25_row"].to_numpy(np.int32),
            "rel_xyz": rel,
            "av": pd.to_numeric(tdf["AV"], errors="coerce").to_numpy(float),
            "original_centre": centre,
        })
        selected_rows.extend(idx.tolist())

    if sum(len(t["indices"]) for t in templates) != 697:
        raise RuntimeError("Template total is not 697; wrong inputs?")
    field_df = df.drop(index=np.unique(selected_rows)).reset_index(drop=True)
    if len(field_df) != 24706 - 697:
        raise RuntimeError("Field row count mismatch")
    return LoadedData(df, field_df, templates, labels)


def make_nominal_field(data: LoadedData, field_id: int) -> dict:
    """Reproduce exactly the outer empirical-template field geometry."""
    rng = np.random.default_rng(OUTER_SEED + field_id)
    field_df = data.field_df
    templates = data.templates

    field_xyz = field_df[["X", "Y", "Z"]].to_numpy(float)
    field_av = pd.to_numeric(field_df["AV"], errors="coerce").to_numpy(float)
    field_q25_rows = field_df["q25_row"].to_numpy(np.int32)
    field_dist = np.linalg.norm(field_xyz, axis=1)
    anchor_pool = sample_anchor_indices(field_xyz, field_dist, rng)

    rotations = [Rotation.random(random_state=rng).as_matrix() for _ in templates]
    centres: list[np.ndarray] = []
    placed: list[np.ndarray] = []
    placement_rows: list[dict] = []

    overlap_direction = field_xyz[int(rng.choice(anchor_pool))]
    overlap_direction = overlap_direction / np.linalg.norm(overlap_direction)
    d0 = float(rng.uniform(600.0, 820.0))
    delta = float(rng.uniform(CFG["overlap_distance_min_pc"], CFG["overlap_distance_max_pc"]))
    overlap_dists = [d0, d0 + delta]

    for idx, template in enumerate(templates):
        rel = template["rel_xyz"]
        rot = rotations[idx]
        if idx < 2:
            centre = overlap_direction * overlap_dists[idx]
            xyz = place_template(rel, centre, rot)
            tries = 0
            while (
                np.linalg.norm(xyz, axis=1).max() > 1000.0
                or np.linalg.norm(xyz, axis=1).min() < 30.0
            ) and tries < 100:
                rot = Rotation.random(random_state=rng).as_matrix()
                xyz = place_template(rel, centre, rot)
                tries += 1
        else:
            xyz = None
            for _ in range(1000):
                centre = field_xyz[int(rng.choice(anchor_pool))].copy()
                if all(np.linalg.norm(centre - c) >= CFG["centre_min_separation_pc"] for c in centres):
                    rot = Rotation.random(random_state=rng).as_matrix()
                    trial = place_template(rel, centre, rot)
                    radii = np.linalg.norm(trial, axis=1)
                    if radii.min() >= 30.0 and radii.max() <= 1000.0:
                        xyz = trial
                        break
            if xyz is None:
                raise RuntimeError(f"Could not place template {idx}")
        centres.append(centre)
        placed.append(xyz)
        placement_rows.append({
            "field_id": field_id,
            "template_label": template["label"],
            "template_gid": template["gid"],
            "n_stars": len(xyz),
            "is_overlap_member": idx < 2,
            "centre_x_pc": float(centre[0]),
            "centre_y_pc": float(centre[1]),
            "centre_z_pc": float(centre[2]),
            "centre_distance_pc": float(np.linalg.norm(centre)),
            "rotation_matrix": json.dumps(rot.tolist()),
            "median_av": float(np.nanmedian(template["av"])),
        })

    cluster_xyz = np.vstack(placed)
    cluster_av = np.concatenate([t["av"] for t in templates])
    cluster_truth = np.concatenate([np.repeat(t["label"], len(t["av"])) for t in templates])
    cluster_q25_rows = np.concatenate([t["q25_rows"] for t in templates]).astype(np.int32)

    xyz = np.vstack([field_xyz, cluster_xyz])
    av = np.concatenate([field_av, cluster_av])
    truth = np.concatenate([np.repeat("background", len(field_df)), cluster_truth])
    q25_rows = np.concatenate([field_q25_rows, cluster_q25_rows])

    order = rng.permutation(len(xyz))
    xyz = xyz[order]
    av = av[order]
    truth = truth[order]
    q25_rows = q25_rows[order]

    # Permanent synthetic-star IDs: source-Q25 row plus occurrence index is not needed,
    # because each source row appears exactly once in a nominal field.
    star_id = q25_rows.copy()
    if len(np.unique(star_id)) != len(star_id):
        raise RuntimeError("Synthetic star IDs are not unique")

    return {
        "xyz": xyz,
        "av": av,
        "truth": truth,
        "q25_rows": q25_rows,
        "star_id": star_id,
        "placements": pd.DataFrame(placement_rows),
        "overlap_pair": [templates[0]["label"], templates[1]["label"]],
        "overlap_distance_separation_pc": delta,
    }


def synthetic_percentiles_from_q25_fractional_uncertainty(
    full_df: pd.DataFrame,
    q25_rows: np.ndarray,
    nominal_xyz: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Transfer each source star's Q25 fractional asymmetric uncertainty to its MC distance."""
    src = full_df.iloc[q25_rows]
    source_d16 = pd.to_numeric(src["b_Dist"], errors="coerce").to_numpy(float)
    source_d50 = pd.to_numeric(src["Dist"], errors="coerce").to_numpy(float)
    source_d84 = pd.to_numeric(src["B_Dist"], errors="coerce").to_numpy(float)
    lower_fraction = (source_d50 - source_d16) / source_d50
    upper_fraction = (source_d84 - source_d50) / source_d50

    d50 = np.linalg.norm(nominal_xyz, axis=1)
    d16 = d50 * (1.0 - lower_fraction)
    d84 = d50 * (1.0 + upper_fraction)
    if np.any(d16 <= 0) or np.any(d84 < d50):
        raise RuntimeError("Invalid synthetic percentile construction")
    return d16, d50, d84


def sample_percentile_matched_distances(
    rng: np.random.Generator,
    d16: np.ndarray,
    d50: np.ndarray,
    d84: np.ndarray,
) -> np.ndarray:
    """Two-sided normal with exact 16th/50th/84th percentiles."""
    z = rng.standard_normal(d50.size)
    q84 = norm.ppf(0.84)
    q16_abs = abs(norm.ppf(0.16))
    lower_scale = (d50 - d16) / q16_abs
    upper_scale = (d84 - d50) / q84
    draws = np.where(z < 0.0, d50 + z * lower_scale, d50 + z * upper_scale)
    invalid = draws <= 0.0
    attempts = 0
    while np.any(invalid):
        attempts += 1
        if attempts > 100:
            raise RuntimeError("Could not obtain positive distance draws")
        z_new = rng.standard_normal(int(invalid.sum()))
        lo = lower_scale[invalid]
        hi = upper_scale[invalid]
        med = d50[invalid]
        draws[invalid] = np.where(z_new < 0.0, med + z_new * lo, med + z_new * hi)
        invalid = draws <= 0.0
    return draws


def evaluate_truth_groups(
    groups: list[np.ndarray],
    truth: np.ndarray,
    template_labels: list[str],
    overlap_pair: list[str],
) -> tuple[pd.DataFrame, dict]:
    group_sizes = np.array([len(g) for g in groups], dtype=int)
    rows = []
    best_ids: dict[str, int] = {}
    for label in template_labels:
        true_members = np.flatnonzero(truth == label)
        true_size = len(true_members)
        best_j = 0.0
        best = -1
        best_inter = 0
        # Match the authoritative MC1000 recovery rule: select the recovered
        # subtree containing the greatest number of true injected members.
        # np.argmax in the production code resolves ties by first occurrence.
        for gid, group in enumerate(groups):
            inter = int(np.intersect1d(true_members, group, assume_unique=True).size)
            if inter > best_inter:
                best, best_inter = gid, inter
        if best >= 0:
            union = true_size + len(groups[best]) - best_inter
            best_j = best_inter / union if union else 0.0
        if best < 0:
            recovered = 0
            gsize = 0
            jaccard = 0.0
            completeness = 0.0
            purity = 0.0
        else:
            recovered = best_inter
            gsize = int(group_sizes[best])
            jaccard = float(best_j)
            completeness = recovered / true_size
            purity = recovered / gsize
        detected = recovered >= max(CFG["min_group_size"], math.ceil(0.5 * true_size))
        best_ids[label] = best
        rows.append({
            "template_label": label,
            "true_size": true_size,
            "best_group_id": best,
            "recovered_members": recovered,
            "best_group_size": gsize,
            "jaccard": jaccard,
            "completeness": completeness,
            "purity": purity,
            "detected_ge50pct": bool(detected),
        })
    tab = pd.DataFrame(rows)
    a, b = overlap_pair
    overlap_resolved = (
        bool(tab.loc[tab.template_label == a, "detected_ge50pct"].iloc[0])
        and bool(tab.loc[tab.template_label == b, "detected_ge50pct"].iloc[0])
        and best_ids[a] >= 0
        and best_ids[b] >= 0
        and best_ids[a] != best_ids[b]
    )
    summary = {
        "templates_detected": int(tab["detected_ge50pct"].sum()),
        "mean_completeness": float(tab["completeness"].mean()),
        "mean_purity": float(tab["purity"].mean()),
        "mean_jaccard": float(tab["jaccard"].mean()),
        "overlap_pair_resolved": bool(overlap_resolved),
    }
    return tab, summary


def match_nominal_groups(
    nominal_groups: list[np.ndarray],
    perturbed_groups: list[np.ndarray],
    perturbed_group_id: np.ndarray,
) -> tuple[pd.DataFrame, list[np.ndarray]]:
    """Independently match every nominal retained group by maximum Jaccard."""
    p_sizes = np.array([len(g) for g in perturbed_groups], dtype=int)
    rows = []
    hit_vectors: list[np.ndarray] = []

    for nominal_gid, members in enumerate(nominal_groups):
        assigned = perturbed_group_id[members]
        candidate_ids = np.unique(assigned[assigned >= 0])
        best_gid = -1
        best_j = -1.0
        best_inter = 0
        for pgid in candidate_ids:
            inter = int(np.sum(assigned == pgid))
            union = len(members) + int(p_sizes[pgid]) - inter
            j = inter / union if union else 0.0
            if j > best_j:
                best_j = j
                best_gid = int(pgid)
                best_inter = inter
        if best_gid < 0:
            best_j = 0.0
            best_inter = 0
            best_size = 0
            hits = np.zeros(len(members), dtype=np.int8)
        else:
            best_size = int(p_sizes[best_gid])
            hits = (assigned == best_gid).astype(np.int8)

        completeness = best_inter / len(members)
        purity = best_inter / best_size if best_size else 0.0
        threshold = max(3, math.ceil(0.10 * len(members)))
        fragment_count = int(sum(np.sum(assigned == pgid) >= threshold for pgid in candidate_ids))
        rows.append({
            "nominal_group_id": nominal_gid,
            "nominal_size": len(members),
            "best_perturbed_group_id": best_gid,
            "best_jaccard": float(best_j),
            "completeness": float(completeness),
            "purity": float(purity),
            "fragment_count": fragment_count,
            "recovered_j_ge_0_50": bool(best_j >= 0.50),
            "recovered_j_ge_0_75": bool(best_j >= 0.75),
            "recovered_j_ge_0_90": bool(best_j >= 0.90),
        })
        hit_vectors.append(hits)
    return pd.DataFrame(rows), hit_vectors


def field_complete(field_dir: Path, expected_inner: int) -> bool:
    marker = field_dir / "field_summary.json"
    if not marker.exists():
        return False
    try:
        data = json.loads(marker.read_text(encoding="utf-8"))
        return int(data.get("uncertainty_realisations", -1)) == expected_inner
    except Exception:
        return False


def process_field(
    field_id: int,
    input_dir_str: str,
    output_dir_str: str,
    uncertainty_runs: int,
    min_group_size: int,
    overwrite: bool,
) -> dict:
    input_dir = Path(input_dir_str)
    output_dir = Path(output_dir_str)
    field_dir = output_dir / "fields" / f"field_{field_id:04d}"
    field_dir.mkdir(parents=True, exist_ok=True)
    if field_complete(field_dir, uncertainty_runs) and not overwrite:
        return json.loads((field_dir / "field_summary.json").read_text(encoding="utf-8"))

    started = time.time()
    data = load_data(input_dir)
    nominal = make_nominal_field(data, field_id)
    xyz0 = nominal["xyz"]
    radii = np.linalg.norm(xyz0, axis=1)
    directions = xyz0 / radii[:, None]
    d16, d50, d84 = synthetic_percentiles_from_q25_fractional_uncertainty(
        data.full_df, nominal["q25_rows"], xyz0
    )

    baseline = run_xmst(xyz0, min_group_size)
    nominal_groups = baseline["groups"]
    baseline_truth, baseline_truth_summary = evaluate_truth_groups(
        nominal_groups, nominal["truth"], data.template_labels, nominal["overlap_pair"]
    )
    baseline_truth.insert(0, "field_id", field_id)
    baseline_truth.to_csv(field_dir / "baseline_template_recovery.csv", index=False)
    nominal["placements"].to_csv(field_dir / "placements.csv", index=False)

    n_nominal = len(nominal_groups)
    jaccard_matrix = np.empty((uncertainty_runs, n_nominal), dtype=np.float32)
    comp_matrix = np.empty((uncertainty_runs, n_nominal), dtype=np.float32)
    purity_matrix = np.empty((uncertainty_runs, n_nominal), dtype=np.float32)
    rec50_matrix = np.empty((uncertainty_runs, n_nominal), dtype=bool)
    rec75_matrix = np.empty((uncertainty_runs, n_nominal), dtype=bool)
    rec90_matrix = np.empty((uncertainty_runs, n_nominal), dtype=bool)
    fragment_matrix = np.empty((uncertainty_runs, n_nominal), dtype=np.int16)
    membership_hits = [np.zeros(len(g), dtype=np.int32) for g in nominal_groups]

    run_rows: list[dict] = []
    template_rows: list[pd.DataFrame] = []

    for inner in range(uncertainty_runs):
        ss = np.random.SeedSequence([INNER_SEED, field_id, inner])
        rng = np.random.default_rng(ss)
        sampled = sample_percentile_matched_distances(rng, d16, d50, d84)
        perturbed_xyz = directions * sampled[:, None]
        pert = run_xmst(perturbed_xyz, min_group_size)

        matches, hit_vectors = match_nominal_groups(nominal_groups, pert["groups"], pert["group_id"])
        jaccard_matrix[inner, :] = matches["best_jaccard"].to_numpy(float)
        comp_matrix[inner, :] = matches["completeness"].to_numpy(float)
        purity_matrix[inner, :] = matches["purity"].to_numpy(float)
        rec50_matrix[inner, :] = matches["recovered_j_ge_0_50"].to_numpy(bool)
        rec75_matrix[inner, :] = matches["recovered_j_ge_0_75"].to_numpy(bool)
        rec90_matrix[inner, :] = matches["recovered_j_ge_0_90"].to_numpy(bool)
        fragment_matrix[inner, :] = matches["fragment_count"].to_numpy(int)
        for gid, hits in enumerate(hit_vectors):
            membership_hits[gid] += hits

        truth_tab, truth_summary = evaluate_truth_groups(
            pert["groups"], nominal["truth"], data.template_labels, nominal["overlap_pair"]
        )
        truth_tab.insert(0, "uncertainty_run", inner)
        truth_tab.insert(0, "field_id", field_id)
        template_rows.append(truth_tab)

        run_rows.append({
            "field_id": field_id,
            "uncertainty_run": inner,
            "nominal_jenks_scale_pc": baseline["jenks_scale_pc"],
            "jenks_scale_pc": pert["jenks_scale_pc"],
            "delta_jenks_pc": pert["jenks_scale_pc"] - baseline["jenks_scale_pc"],
            "percolation_cut_pc": pert["percolation_cut_pc"],
            "n_retained_groups": len(pert["groups"]),
            "n_grouped_stars": int(np.sum(pert["group_id"] >= 0)),
            "mean_nominal_group_jaccard": float(matches["best_jaccard"].mean()),
            "median_nominal_group_jaccard": float(matches["best_jaccard"].median()),
            "fraction_nominal_groups_j_ge_0_50": float(matches["recovered_j_ge_0_50"].mean()),
            "fraction_nominal_groups_j_ge_0_75": float(matches["recovered_j_ge_0_75"].mean()),
            "fraction_nominal_groups_j_ge_0_90": float(matches["recovered_j_ge_0_90"].mean()),
            "mean_nominal_group_completeness": float(matches["completeness"].mean()),
            "mean_nominal_group_purity": float(matches["purity"].mean()),
            "templates_detected": truth_summary["templates_detected"],
            "template_mean_completeness": truth_summary["mean_completeness"],
            "template_mean_purity": truth_summary["mean_purity"],
            "template_mean_jaccard": truth_summary["mean_jaccard"],
            "overlap_pair_resolved": truth_summary["overlap_pair_resolved"],
            "distance_ratio_median": float(np.median(sampled / d50)),
            "distance_ratio_p16": float(np.quantile(sampled / d50, 0.16)),
            "distance_ratio_p84": float(np.quantile(sampled / d50, 0.84)),
        })

    runs_df = pd.DataFrame(run_rows)
    runs_df.to_csv(field_dir / "uncertainty_runs.csv", index=False)
    if template_rows:
        pd.concat(template_rows, ignore_index=True).to_csv(
            field_dir / "template_recovery_all_uncertainty_runs.csv", index=False
        )

    # Aggregate group stability for this nominal field.
    group_rows = []
    for gid, members in enumerate(nominal_groups):
        j = jaccard_matrix[:, gid]
        c = comp_matrix[:, gid]
        p = purity_matrix[:, gid]
        group_rows.append({
            "field_id": field_id,
            "nominal_group_id": gid,
            "nominal_size": len(members),
            "mean_best_jaccard": float(np.mean(j)),
            "median_best_jaccard": float(np.median(j)),
            "jaccard_p16": float(np.quantile(j, 0.16)),
            "jaccard_p84": float(np.quantile(j, 0.84)),
            "recovery_fraction_j_ge_0_50": float(np.mean(rec50_matrix[:, gid])),
            "recovery_fraction_j_ge_0_75": float(np.mean(rec75_matrix[:, gid])),
            "recovery_fraction_j_ge_0_90": float(np.mean(rec90_matrix[:, gid])),
            "median_completeness": float(np.median(c)),
            "median_purity": float(np.median(p)),
            "mean_fragment_count": float(np.mean(fragment_matrix[:, gid])),
        })
    pd.DataFrame(group_rows).to_csv(field_dir / "nominal_group_stability.csv", index=False)

    # Aggregate star persistence for stars that belong to nominal retained groups.
    star_rows = []
    for gid, members in enumerate(nominal_groups):
        scores = membership_hits[gid] / max(1, uncertainty_runs)
        for node, score in zip(members, scores):
            star_rows.append({
                "field_id": field_id,
                "nominal_group_id": gid,
                "synthetic_node": int(node),
                "source_q25_row": int(nominal["q25_rows"][node]),
                "true_component": str(nominal["truth"][node]),
                "membership_persistence": float(score),
            })
    pd.DataFrame(star_rows).to_csv(field_dir / "star_membership_persistence.csv", index=False)

    summary = {
        "field_id": field_id,
        "outer_seed": OUTER_SEED + field_id,
        "inner_master_seed": INNER_SEED,
        "uncertainty_realisations": uncertainty_runs,
        "minimum_group_size": min_group_size,
        "n_stars": len(xyz0),
        "n_nominal_groups": n_nominal,
        "n_nominal_grouped_stars": int(np.sum(baseline["group_id"] >= 0)),
        "nominal_percolation_cut_pc": baseline["percolation_cut_pc"],
        "nominal_jenks_scale_pc": baseline["jenks_scale_pc"],
        "nominal_templates_detected": baseline_truth_summary["templates_detected"],
        "nominal_template_mean_completeness": baseline_truth_summary["mean_completeness"],
        "nominal_template_mean_purity": baseline_truth_summary["mean_purity"],
        "nominal_template_mean_jaccard": baseline_truth_summary["mean_jaccard"],
        "nominal_overlap_pair_resolved": baseline_truth_summary["overlap_pair_resolved"],
        "mc_jenks_mean_pc": float(runs_df["jenks_scale_pc"].mean()),
        "mc_jenks_median_pc": float(runs_df["jenks_scale_pc"].median()),
        "mc_jenks_sd_pc": float(runs_df["jenks_scale_pc"].std(ddof=1)) if uncertainty_runs > 1 else 0.0,
        "mc_delta_jenks_mean_pc": float(runs_df["delta_jenks_pc"].mean()),
        "mc_delta_jenks_abs_median_pc": float(runs_df["delta_jenks_pc"].abs().median()),
        "mean_of_run_mean_jaccard": float(runs_df["mean_nominal_group_jaccard"].mean()),
        "mean_fraction_nominal_groups_j_ge_0_50": float(runs_df["fraction_nominal_groups_j_ge_0_50"].mean()),
        "mean_fraction_nominal_groups_j_ge_0_75": float(runs_df["fraction_nominal_groups_j_ge_0_75"].mean()),
        "mean_fraction_nominal_groups_j_ge_0_90": float(runs_df["fraction_nominal_groups_j_ge_0_90"].mean()),
        "mean_template_completeness_under_uncertainty": float(runs_df["template_mean_completeness"].mean()),
        "mean_template_purity_under_uncertainty": float(runs_df["template_mean_purity"].mean()),
        "mean_template_jaccard_under_uncertainty": float(runs_df["template_mean_jaccard"].mean()),
        "fraction_runs_all_8_templates_detected": float(np.mean(runs_df["templates_detected"] == 8)),
        "fraction_runs_overlap_pair_resolved": float(runs_df["overlap_pair_resolved"].mean()),
        "elapsed_seconds": time.time() - started,
    }
    (field_dir / "field_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    return summary


def aggregate(output_dir: Path, field_ids: Iterable[int], uncertainty_runs: int) -> None:
    summaries = []
    run_frames = []
    group_frames = []
    star_frames = []
    template_frames = []
    for field_id in field_ids:
        d = output_dir / "fields" / f"field_{field_id:04d}"
        if not field_complete(d, uncertainty_runs):
            continue
        summaries.append(json.loads((d / "field_summary.json").read_text(encoding="utf-8")))
        run_frames.append(pd.read_csv(d / "uncertainty_runs.csv"))
        group_frames.append(pd.read_csv(d / "nominal_group_stability.csv"))
        star_frames.append(pd.read_csv(d / "star_membership_persistence.csv"))
        template_frames.append(pd.read_csv(d / "template_recovery_all_uncertainty_runs.csv"))

    if not summaries:
        return
    summary_df = pd.DataFrame(summaries).sort_values("field_id")
    runs_df = pd.concat(run_frames, ignore_index=True)
    groups_df = pd.concat(group_frames, ignore_index=True)
    stars_df = pd.concat(star_frames, ignore_index=True)
    templates_df = pd.concat(template_frames, ignore_index=True)

    summary_df.to_csv(output_dir / "outer_field_summaries.csv", index=False)
    runs_df.to_csv(output_dir / "all_uncertainty_runs.csv", index=False)
    groups_df.to_csv(output_dir / "all_nominal_group_stability.csv", index=False)
    stars_df.to_csv(output_dir / "all_star_membership_persistence.csv", index=False)
    templates_df.to_csv(output_dir / "all_template_recovery.csv", index=False)

    # Baseline reproduction quantities over completed outer fields.
    baseline = {
        "completed_outer_fields": int(len(summary_df)),
        "uncertainty_realisations_per_field": int(uncertainty_runs),
        "total_uncertainty_xmst_runs": int(len(runs_df)),
        "baseline_median_jenks_pc": float(summary_df["nominal_jenks_scale_pc"].median()),
        "baseline_sd_jenks_pc": float(summary_df["nominal_jenks_scale_pc"].std(ddof=1)) if len(summary_df) > 1 else 0.0,
        "baseline_mean_template_completeness": float(summary_df["nominal_template_mean_completeness"].mean()),
        "baseline_mean_template_purity": float(summary_df["nominal_template_mean_purity"].mean()),
        "perturbed_jenks_median_pc": float(runs_df["jenks_scale_pc"].median()),
        "perturbed_jenks_sd_pc": float(runs_df["jenks_scale_pc"].std(ddof=1)) if len(runs_df) > 1 else 0.0,
        "delta_jenks_median_pc": float(runs_df["delta_jenks_pc"].median()),
        "absolute_delta_jenks_median_pc": float(runs_df["delta_jenks_pc"].abs().median()),
        "mean_nominal_group_jaccard": float(runs_df["mean_nominal_group_jaccard"].mean()),
        "median_nominal_group_jaccard": float(groups_df["median_best_jaccard"].median()),
        "mean_group_recovery_fraction_j_ge_0_50": float(groups_df["recovery_fraction_j_ge_0_50"].mean()),
        "mean_group_recovery_fraction_j_ge_0_75": float(groups_df["recovery_fraction_j_ge_0_75"].mean()),
        "mean_group_recovery_fraction_j_ge_0_90": float(groups_df["recovery_fraction_j_ge_0_90"].mean()),
        "median_star_membership_persistence": float(stars_df["membership_persistence"].median()),
        "star_fraction_persistence_ge_0_50": float(np.mean(stars_df["membership_persistence"] >= 0.50)),
        "star_fraction_persistence_ge_0_80": float(np.mean(stars_df["membership_persistence"] >= 0.80)),
        "star_fraction_persistence_ge_0_90": float(np.mean(stars_df["membership_persistence"] >= 0.90)),
        "mean_template_completeness_under_uncertainty": float(templates_df["completeness"].mean()),
        "mean_template_purity_under_uncertainty": float(templates_df["purity"].mean()),
        "mean_template_jaccard_under_uncertainty": float(templates_df["jaccard"].mean()),
        "fraction_template_cases_detected": float(templates_df["detected_ge50pct"].mean()),
    }
    (output_dir / "overall_summary.json").write_text(json.dumps(baseline, indent=2), encoding="utf-8")

    report = [
        "# XMST nested distance-uncertainty experiment",
        "",
        f"Completed outer fields: **{baseline['completed_outer_fields']}**",
        f"Uncertainty realisations per field: **{uncertainty_runs}**",
        f"Total perturbed XMST reconstructions: **{baseline['total_uncertainty_xmst_runs']}**",
        "",
        "## Baseline reproduction",
        f"- Median nominal Jenks scale: {baseline['baseline_median_jenks_pc']:.6f} pc",
        f"- Between-field nominal Jenks SD: {baseline['baseline_sd_jenks_pc']:.6f} pc",
        f"- Mean nominal template completeness: {baseline['baseline_mean_template_completeness']:.6f}",
        f"- Mean nominal template purity: {baseline['baseline_mean_template_purity']:.6f}",
        "",
        "## Distance-uncertainty result",
        f"- Median perturbed Jenks scale: {baseline['perturbed_jenks_median_pc']:.6f} pc",
        f"- Perturbed Jenks SD: {baseline['perturbed_jenks_sd_pc']:.6f} pc",
        f"- Median Jenks shift relative to each nominal field: {baseline['delta_jenks_median_pc']:.6f} pc",
        f"- Median absolute Jenks shift: {baseline['absolute_delta_jenks_median_pc']:.6f} pc",
        f"- Mean nominal-group Jaccard: {baseline['mean_nominal_group_jaccard']:.6f}",
        f"- Median group median-Jaccard: {baseline['median_nominal_group_jaccard']:.6f}",
        f"- Mean group recovery fraction (J >= 0.50): {baseline['mean_group_recovery_fraction_j_ge_0_50']:.6f}",
        f"- Mean group recovery fraction (J >= 0.75): {baseline['mean_group_recovery_fraction_j_ge_0_75']:.6f}",
        f"- Mean group recovery fraction (J >= 0.90): {baseline['mean_group_recovery_fraction_j_ge_0_90']:.6f}",
        f"- Median star membership persistence: {baseline['median_star_membership_persistence']:.6f}",
        f"- Stars with persistence >= 0.80: {baseline['star_fraction_persistence_ge_0_80']:.6f}",
        f"- Mean injected-template completeness under uncertainty: {baseline['mean_template_completeness_under_uncertainty']:.6f}",
        f"- Mean injected-template purity under uncertainty: {baseline['mean_template_purity_under_uncertainty']:.6f}",
        f"- Fraction of template cases satisfying >=50% detection: {baseline['fraction_template_cases_detected']:.6f}",
        "",
        "## Sampling definition",
        "For every synthetic star, the lower and upper *fractional* Q25 distance uncertainties are preserved.",
        "The moved star's nominal MC radius is treated as its median distance, and synthetic 16th/84th percentiles",
        "are constructed by applying the source star's fractional lower/upper uncertainty. Distances are then drawn",
        "from the same percentile-matched two-sided normal approximation used in the Q25 stability code.",
    ]
    (output_dir / "RESULTS_SUMMARY.md").write_text("\n".join(report) + "\n", encoding="utf-8")


def parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--input-dir", type=Path, default=DEFAULT_INPUT)
    ap.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT)
    ap.add_argument("--field-start", type=int, default=0)
    ap.add_argument("--fields", type=int, default=1, help="Number of nominal outer fields")
    ap.add_argument("--uncertainty-runs", type=int, default=100, help="Distance realisations per nominal field")
    ap.add_argument("--workers", type=int, default=max(1, min(8, (os.cpu_count() or 4) // 2)))
    ap.add_argument("--min-group-size", type=int, default=10)
    ap.add_argument("--overwrite", action="store_true")
    return ap.parse_args()


def main() -> int:
    args = parse_args()
    if args.fields < 1 or args.uncertainty_runs < 1:
        raise ValueError("--fields and --uncertainty-runs must be positive")
    field_ids = list(range(args.field_start, args.field_start + args.fields))
    args.output_dir.mkdir(parents=True, exist_ok=True)
    (args.output_dir / "fields").mkdir(parents=True, exist_ok=True)

    config = {
        "field_ids": field_ids,
        "uncertainty_runs_per_field": args.uncertainty_runs,
        "outer_seed": OUTER_SEED,
        "inner_master_seed": INNER_SEED,
        "minimum_group_size": args.min_group_size,
        "distance_sampler": "percentile-matched two-sided normal",
        "moved_template_uncertainty_rule": "preserve source Q25 fractional lower/upper distance uncertainty",
        "spatial_only": True,
        "template_matching_rule": "maximum injected-member overlap (matches authoritative MC1000 recovery rule)",
        "nominal_group_matching_rule": "maximum Jaccard similarity",
        "input_dir": str(args.input_dir.resolve()),
        "software": {
            "python": sys.version,
            "platform": platform.platform(),
            "numpy": np.__version__,
            "pandas": pd.__version__,
            "scipy": scipy.__version__,
        },
    }
    (args.output_dir / "run_config.json").write_text(json.dumps(config, indent=2), encoding="utf-8")

    workers = max(1, min(args.workers, len(field_ids)))
    print("XMST nested distance-uncertainty experiment", flush=True)
    print(f"Outer fields: {field_ids[0]}..{field_ids[-1]} ({len(field_ids)})", flush=True)
    print(f"Uncertainty realisations per field: {args.uncertainty_runs}", flush=True)
    print(f"Total perturbed XMST runs: {len(field_ids) * args.uncertainty_runs:,}", flush=True)
    print(f"Workers: {workers}", flush=True)
    print("A_V pruning: OFF (spatial stability experiment)", flush=True)

    started = time.time()
    if workers == 1:
        for i, field_id in enumerate(field_ids, 1):
            s = process_field(
                field_id, str(args.input_dir), str(args.output_dir),
                args.uncertainty_runs, args.min_group_size, args.overwrite,
            )
            print(
                f"[{i}/{len(field_ids)}] field {field_id:04d}: "
                f"nominal J={s['nominal_jenks_scale_pc']:.3f} pc, "
                f"perturbed median J={s['mc_jenks_median_pc']:.3f} pc, "
                f"mean group Jaccard={s['mean_of_run_mean_jaccard']:.3f}, "
                f"{s['elapsed_seconds']:.1f}s",
                flush=True,
            )
    else:
        # Prevent BLAS/OpenMP oversubscription in each worker.
        for k in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
            os.environ[k] = "1"
        with ProcessPoolExecutor(max_workers=workers) as ex:
            futures = {
                ex.submit(
                    process_field,
                    field_id, str(args.input_dir), str(args.output_dir),
                    args.uncertainty_runs, args.min_group_size, args.overwrite,
                ): field_id
                for field_id in field_ids
            }
            done = 0
            for fut in as_completed(futures):
                field_id = futures[fut]
                s = fut.result()
                done += 1
                print(
                    f"[{done}/{len(field_ids)}] field {field_id:04d}: "
                    f"nominal J={s['nominal_jenks_scale_pc']:.3f} pc, "
                    f"perturbed median J={s['mc_jenks_median_pc']:.3f} pc, "
                    f"mean group Jaccard={s['mean_of_run_mean_jaccard']:.3f}",
                    flush=True,
                )

    aggregate(args.output_dir, field_ids, args.uncertainty_runs)
    elapsed = time.time() - started
    print(f"Completed in {elapsed / 60.0:.1f} min", flush=True)
    print(f"Summary: {args.output_dir / 'RESULTS_SUMMARY.md'}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
