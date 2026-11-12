#!/usr/bin/env python3
"""
hotspot_terminal_filter.py

Screen predicted binder-target complexes (e.g. VHH/nanobody - GPCR) by
(a) counting how many predefined target "hotspot" residues the binder
    contacts, and
(b) rejecting designs whose binder N- or C-terminus buries into the hotspot
    pocket (terminal tags / Fc fusions would then block binding).

Ranking rule: hotspot count (desc), then total interface contacts (desc).

Input
-----
A directory of PDB files, each containing exactly two chains:
  --chain_vhh : the binder chain (default A)
  --chain_tgt : the target chain  (default B)

Hotspots
--------
Given as a comma list. Each entry is either:
  <resnum>                -> use all heavy atoms of that target residue
  <resnum>:<atom>[,atom]  -> use only the listed atoms (side-chain specific)
Example (NK2R, 15 NKA-activation contacts):
  86:ND2,90:ND2,93:OH,94:CB,96:NE2,114:CG1,CG2,CD1,117:SD,CE,
  175:OD1,OD2,180:NZ,181:SG,269:OH,270:CZ,CE1,CE2,280:OH,289:OH,293:CZ,CE1,CE2

Terminal filter
---------------
A design FAILS if the binder's N-terminal window (residues 1..term_len) or
C-terminal window (last term_len residues) contacts >= term_max_hs distinct
hotspots. Rationale: affinity tags (His-tag, FLAG) or C-terminal Fc fusions
are attached at the termini; if a terminus sits inside the binding pocket the
fusion sterically blocks engagement.

Output
------
TSV with one row per design:
  rank design vhh_len n_hotspots hotspot_list n_contacts
       nterm_hs_count cterm_hs_count terminal_pass fail_reason
Survivors (terminal_pass=PASS) are ranked first.

Dependencies: numpy, scipy, biopython
"""
from __future__ import annotations

import argparse
import glob
import os
from collections import defaultdict

import numpy as np
from Bio.PDB import PDBParser
from scipy.spatial import cKDTree

HEAVY = ("N", "CA", "C", "O", "CB", "CG", "CD", "CE", "CZ",
         "NH1", "NH2", "ND1", "ND2", "NE", "NE1", "NE2",
         "OG", "OG1", "OH", "OD1", "OD2", "OE1", "OE2", "SG", "SD", "NZ")


def parse_hotspots(spec: str) -> dict:
    """'86:ND2,90,93:OH' -> {86: ['ND2'], 90: [], 93: ['OH']}  ([]=all heavy)."""
    out = {}
    if not spec:
        return out
    for tok in spec.split(","):
        tok = tok.strip()
        if not tok:
            continue
        if ":" in tok:
            num, atoms = tok.split(":", 1)
            out[int(num)] = [a.strip() for a in atoms.split("/") if a.strip()]
        else:
            out[int(tok)] = []
    return out


def heavy_coords(chain, resnum=None, atom_names=None):
    coords = []
    residues = [chain[resnum]] if resnum is not None else chain
    for res in residues:
        for atom in res:
            if atom.get_name() in HEAVY:
                if atom_names is None or not atom_names or atom.get_name() in atom_names:
                    coords.append(atom.get_coord())
    return np.array(coords) if coords else np.zeros((0, 3))


def analyze(pdb_path, chain_vhh, chain_tgt, hotspots, cutoff):
    parser = PDBParser(QUIET=True)
    struct = parser.get_structure("s", pdb_path)
    chains = {}
    for model in struct:
        for ch in model:
            chains[ch.get_id()] = ch
        break
    if chain_vhh not in chains or chain_tgt not in chains:
        return None
    vhh, tgt = chains[chain_vhh], chains[chain_tgt]

    # binder heavy atoms + their residue ids
    v_coords, v_res = [], []
    for atom in vhh.get_atoms():
        if atom.get_name() in HEAVY:
            v_coords.append(atom.get_coord())
            v_res.append(atom.get_parent().id[1])
    if not v_coords:
        return None
    v_coords = np.array(v_coords)
    v_res = np.array(v_res)
    tree_vhh = cKDTree(v_coords)

    # per-hotspot: does binder touch it? which binder residues?
    hs_hits, hs_vhh_res = [], defaultdict(set)
    for resnum, atoms in sorted(hotspots.items()):
        try:
            hc = heavy_coords(tgt, resnum, atoms or None)
        except KeyError:
            continue
        if len(hc) == 0:
            continue
        idx_list = tree_vhh.query_ball_point(hc, r=cutoff)
        touched = set()
        for idxs in idx_list:
            for ii in idxs:
                touched.add(int(v_res[ii]))
        if touched:
            hs_hits.append(resnum)
            hs_vhh_res[resnum] = touched

    # total contacted target residues (any heavy-heavy <= cutoff)
    t_coords, t_res = [], []
    for atom in tgt.get_atoms():
        if atom.get_name() in HEAVY:
            t_coords.append(atom.get_coord())
            t_res.append(atom.get_parent().id[1])
    n_contacts = 0
    if t_coords:
        tree_tgt = cKDTree(np.array(t_coords))
        d, i = tree_tgt.query(v_coords, k=1)
        n_contacts = len({int(t_res[k]) for k, dd in enumerate(d) if dd <= cutoff})

    # terminal windows are computed in main from rmin/rmax + hs_vhh_res
    rmin, rmax = int(v_res.min()), int(v_res.max())
    n_len = rmax - rmin + 1
    return dict(hs_hits=hs_hits, hs_vhh_res=hs_vhh_res, n_contacts=n_contacts,
                rmin=rmin, rmax=rmax, n_len=n_len)


def main():
    ap = argparse.ArgumentParser(
        description="Rank binder-target complexes by hotspot coverage and "
                    "reject terminus-buried designs.")
    ap.add_argument("--pdb_dir", required=True, help="dir with 2-chain PDBs")
    ap.add_argument("--chain_vhh", default="A")
    ap.add_argument("--chain_tgt", default="B")
    ap.add_argument("--hotspots", required=True,
                    help="e.g. '86:ND2,90,93:OH,...' (see docstring)")
    ap.add_argument("--cutoff", type=float, default=4.5,
                    help="heavy-atom contact cutoff (A)")
    ap.add_argument("--term_len", type=int, default=25,
                    help="terminal window length (residues)")
    ap.add_argument("--term_max_hs", type=int, default=2,
                    help="fail if a terminal window touches >= this many hotspots")
    ap.add_argument("--out", default="hotspot_terminal_ranking.tsv")
    args = ap.parse_args()

    hotspots = parse_hotspots(args.hotspots)
    files = sorted(glob.glob(os.path.join(args.pdb_dir, "*.pdb")))
    if not files:
        files = sorted(glob.glob(os.path.join(args.pdb_dir, "**", "*.pdb"),
                                 recursive=True))
    print("PDB files: %d | hotspots: %d" % (len(files), len(hotspots)))

    rows = []
    for pf in files:
        name = os.path.basename(pf)[:-4]
        try:
            r = analyze(pf, args.chain_vhh, args.chain_tgt, hotspots, args.cutoff)
        except Exception:
            r = None
        if r is None:
            continue
        # terminal classification
        nterm_hs, cterm_hs = [], []
        n_lo, n_hi = r["rmin"], r["rmin"] + args.term_len - 1
        c_lo, c_hi = r["rmax"] - args.term_len + 1, r["rmax"]
        for resnum in r["hs_hits"]:
            vset = r["hs_vhh_res"][resnum]
            if any(n_lo <= x <= n_hi for x in vset):
                nterm_hs.append(resnum)
            if any(c_lo <= x <= c_hi for x in vset):
                cterm_hs.append(resnum)
        fail = []
        if len(nterm_hs) >= args.term_max_hs:
            fail.append("Nterm(%d hs)" % len(nterm_hs))
        if len(cterm_hs) >= args.term_max_hs:
            fail.append("Cterm(%d hs)" % len(cterm_hs))
        rows.append(dict(name=name, n_hs=len(r["hs_hits"]),
                         hs_list=r["hs_hits"], n_contacts=r["n_contacts"],
                         vhh_len=r["n_len"], n_nterm=len(nterm_hs),
                         n_cterm=len(cterm_hs), passed=not fail,
                         fail=";".join(fail)))

    rows.sort(key=lambda x: (not x["passed"], -x["n_hs"], -x["n_contacts"]))

    with open(args.out, "w") as f:
        f.write("rank\tname\tvhh_len\tn_hotspots\thotspot_list\tn_contacts\t"
                "nterm_hs\tcterm_hs\tterminal_pass\tfail_reason\n")
        for i, r in enumerate(rows, 1):
            f.write("%d\t%s\t%d\t%d\t%s\t%d\t%d\t%d\t%s\t%s\n" % (
                i, r["name"], r["vhh_len"], r["n_hs"],
                ",".join(map(str, r["hs_list"])), r["n_contacts"],
                r["n_nterm"], r["n_cterm"],
                "PASS" if r["passed"] else "FAIL", r["fail"]))

    n_pass = sum(1 for r in rows if r["passed"])
    print("Analyzed %d | PASS %d | FAIL %d" % (len(rows), n_pass, len(rows) - n_pass))
    print("\nTop 20 survivors:")
    print("%-4s %-20s %-6s %-8s %s" % ("Rnk", "Design", "Len", "Hs", "Hotspots"))
    k = 0
    for r in rows:
        if not r["passed"]:
            continue
        k += 1
        if k > 20:
            break
        print("%-4d %-20s %-6d %-8s %s" % (k, r["name"], r["vhh_len"],
                                           "%d/%d" % (r["n_hs"], len(hotspots)),
                                           ",".join(map(str, r["hs_list"]))))
    print("\nSaved: %s" % args.out)


if __name__ == "__main__":
    main()
