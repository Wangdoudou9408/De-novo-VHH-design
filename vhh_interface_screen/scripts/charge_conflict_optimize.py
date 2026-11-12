#!/usr/bin/env python3
"""
charge_conflict_optimize.py

Detect same-charge electrostatic repulsions at a binder-target protein
interface and propose single-point mutations on the BINDER side only that
relieve each conflict, for downstream local sequence optimization.
(The target is treated as fixed; this tool suggests, it does not filter.)

Method
------
1. Charged groups are defined at atom level:
       +  Lys(NZ), Arg(NE/NH1/NH2), His(ND1/NE2, weak +0.5)
       -  Asp(OD1/OD2), Glu(OE1/OE2)
2. Interface binder residues = any binder heavy atom within --iface_cutoff
   of a target heavy atom (complete set via KD-tree query_ball_point).
3. For every binder/target charged-group pair:
       opposite charge, d <= --saltbridge      -> favourable salt bridge (kept)
       same charge,     d <= --conflict_strong -> STRONG repulsion
       same charge,     d <= --conflict_mid    -> MID repulsion
4. For each repulsion, candidate binder mutations are scored:
       charge-flip  (K/R->D/E or D/E->K/R)  : turns repulsion into a salt
                                              bridge (highest gain)
       neutralize   (->Q/N/S/A)              : removes the charge (safe)
   Score = gain
         - steric penalty      (side-chain volume mismatch > 60 A^3)
         - hydrophobic penalty (KYTTE-DELFOOR-like scale; disfavour
                                hydrophobic at a charged interface position)
         - buried-charge penalty (new charge at a buried, low-SASA site)
         - protection penalty  (residue contacts a hotspot: flip forbidden,
                                neutralize down-weighted; see --allow-*)
         - CDR penalty         (with --use_anarci: CDR positions never flip)
         - internal-context penalties:
             * losing an existing binder-internal salt bridge
             * creating a new binder-internal same-charge clash
             * creating a new binder-internal salt bridge (CDR risk)
         - new-interface-conflict penalty:
             the proposed charge is re-scanned against ALL nearby target
             charged groups; new same-charge clashes are penalised
5. Hotspot-contacting binder residues are PROTECTED: by default reported but
   skipped; with --allow_protected_neutralize they may only be neutralized
   (never flipped).

Output
------
TSV, one row per charged pair:
  design vhh_res vhh_aa tgt_res tgt_aa pair_type severity dist_A
  suggested_mut strategy score note
`note` is one of: mutate | PROTECTED-neutralize-only | CDR-neutralize-only |
PROTECTED(hotspot-contact): skip | no-safe-mutation | keep (salt bridge).

Dependencies: numpy, scipy, biopython. Optional: ANARCI (CDR/FR split).
"""
import os, argparse, glob
import numpy as np
from Bio.PDB import PDBParser, ShrakeRupley
from scipy.spatial import cKDTree

POS_ATOMS = {"LYS": ["NZ"], "ARG": ["NE", "NH1", "NH2"], "HIS": ["ND1", "NE2"]}
NEG_ATOMS = {"ASP": ["OD1", "OD2"], "GLU": ["OE1", "OE2"]}
ONE = {"LYS": "K", "ARG": "R", "HIS": "H", "ASP": "D", "GLU": "E"}
AA3TO1 = {"ALA":"A","ARG":"R","ASN":"N","ASP":"D","CYS":"C","GLN":"Q","GLU":"E","GLY":"G",
          "HIS":"H","ILE":"I","LEU":"L","LYS":"K","MET":"M","PHE":"F","PRO":"P","SER":"S",
          "THR":"T","TRP":"W","TYR":"Y","VAL":"V"}

VOL = {"A": 88.6, "R": 173.4, "N": 114.1, "D": 111.1, "C": 108.5, "Q": 143.8,
       "E": 138.4, "G": 60.1, "H": 153.2, "I": 166.7, "L": 166.7, "K": 168.6,
       "M": 162.9, "F": 189.9, "P": 112.7, "S": 89.0, "T": 116.1, "W": 227.8,
       "Y": 193.6, "V": 140.0}
HYDRO = {"A": 1.8, "R": -4.5, "N": -3.5, "D": -3.5, "C": 2.5, "Q": -3.5,
         "E": -3.5, "G": -0.4, "H": -3.2, "I": 4.5, "L": 3.8, "K": -3.9,
         "M": 1.9, "F": 2.8, "P": -1.6, "S": -0.8, "T": -0.7, "W": -0.9,
         "Y": -1.3, "V": 4.2}

try:
    from anarci import anarci as _anarci
    HAVE_ANARCI = True
except Exception:
    HAVE_ANARCI = False


def cdr_set_anarci(seq):
    """Return set of 1-based SEQUENCE positions (non-gap index) in CDR (IMGT for VHH)."""
    if not HAVE_ANARCI:
        return set()
    try:
        results = _anarci([("vhh", seq)], scheme="imgt", output=False)
        if not (results and results[0] and results[0][0]):
            return set()
        cdr = set()
        seq_idx = 0                       # actual sequence position, 1-based
        for (pos, aa) in results[0][0][0]:
            if aa == "-":
                continue
            seq_idx += 1                  # only increment on non-gap
            num = pos[0]
            if (27 <= num <= 38) or (56 <= num <= 65) or (105 <= num <= 117):
                cdr.add(seq_idx)          # store sequence index, NOT IMGT number
        return cdr
    except Exception as e:
        print("[warn] ANARCI failed:", e)   # do not swallow silently
        return set()


def charged_groups(chain):
    groups = []
    for res in chain:
        rn = res.get_resname().strip()
        if rn in POS_ATOMS:
            coords = [res[a].get_coord() for a in POS_ATOMS[rn] if a in res]
            if coords:
                groups.append({"res": res.id[1], "aa": ONE[rn],
                               "q": (+0.5 if rn == "HIS" else +1),
                               "coords": np.array(coords)})
        elif rn in NEG_ATOMS:
            coords = [res[a].get_coord() for a in NEG_ATOMS[rn] if a in res]
            if coords:
                groups.append({"res": res.id[1], "aa": ONE[rn], "q": -1,
                               "coords": np.array(coords)})
    return groups


def gdist(c1, c2):
    d = np.sqrt(((c1[:, None, :] - c2[None, :, :]) ** 2).sum(-1))
    return d.min()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pdb_dir", required=True)
    ap.add_argument("--chain_vhh", default="A")
    ap.add_argument("--chain_tgt", default="B")
    ap.add_argument("--hotspots", default="")
    ap.add_argument("--iface_cutoff", type=float, default=5.0)
    ap.add_argument("--conflict_strong", type=float, default=4.0)
    ap.add_argument("--conflict_mid", type=float, default=5.0)
    ap.add_argument("--saltbridge", type=float, default=4.0)
    ap.add_argument("--internal_cutoff", type=float, default=6.0)
    ap.add_argument("--allow_protected_neutralize", action="store_true")
    ap.add_argument("--use_anarci", action="store_true")
    ap.add_argument("--out", default="charge_optimization.tsv")
    args = ap.parse_args()

    hotspots = set(int(x) for x in args.hotspots.split(",") if x.strip())
    files = sorted(glob.glob(os.path.join(args.pdb_dir, "*.pdb")))
    if not files:
        files = sorted(glob.glob(os.path.join(args.pdb_dir, "**", "*.pdb"), recursive=True))
    print("Structures: %d | ANARCI: %s" % (len(files), (HAVE_ANARCI and args.use_anarci)))

    rows = []
    for pf in files:
        name = os.path.basename(pf).replace(".pdb", "")
        try:
            s = PDBParser(QUIET=True).get_structure(name, pf)
        except Exception:
            continue
        ch = {}
        for m in s:
            for c in m:
                ch[c.get_id()] = c
            break
        if args.chain_vhh not in ch or args.chain_tgt not in ch:
            continue
        vhh, tgt = ch[args.chain_vhh], ch[args.chain_tgt]

        try:
            ShrakeRupley().compute(s, level="R")
        except Exception:
            pass

        v_coords, v_res = [], []
        for a in vhh.get_atoms():
            if a.element != "H":
                v_coords.append(a.coord); v_res.append(a.get_parent().id[1])
        t_coords = []
        for a in tgt.get_atoms():
            if a.element != "H":
                t_coords.append(a.coord)
        if not v_coords or not t_coords:
            continue
        v_coords = np.array(v_coords); v_res = np.array(v_res)
        t_coords = np.array(t_coords)
        tree_t = cKDTree(t_coords)
        tree_v = cKDTree(v_coords)

        # FIX 3: complete interface set via query_ball_point
        ball = tree_t.query_ball_point(v_coords, r=args.iface_cutoff)
        iface_vhh_res = set(int(v_res[k]) for k, idxs in enumerate(ball) if idxs)

        # FIX 2: complete hotspot-contact protection via query_ball_point
        hs_contact_vhh = set()
        for h in hotspots:
            try:
                hres = tgt[h]
            except KeyError:
                continue
            hc = np.array([a.get_coord() for a in hres.get_atoms() if a.element != "H"])
            if len(hc) == 0:
                continue
            idx_list = tree_v.query_ball_point(hc, r=args.iface_cutoff)
            for idxs in idx_list:
                for ii in idxs:
                    hs_contact_vhh.add(int(v_res[ii]))

        # CDR mapping (optional)
        vhh_seq = ""
        resnum_to_seqidx = {}
        for r in vhh:
            vhh_seq += AA3TO1.get(r.get_resname().strip(), "X")
            resnum_to_seqidx[r.id[1]] = len(vhh_seq)
        cdr_seqidx = cdr_set_anarci(vhh_seq) if (args.use_anarci and HAVE_ANARCI) else set()
        seqidx_to_res = {v: k for k, v in resnum_to_seqidx.items()}
        cdr_res = set(seqidx_to_res[i] for i in cdr_seqidx if i in seqidx_to_res)

        gv = charged_groups(vhh)
        gt = charged_groups(tgt)

        for g in gv:
            vr = g["res"]
            if vr not in iface_vhh_res:
                continue
            vq, v1 = g["q"], g["aa"]
            for h in gt:
                d = gdist(g["coords"], h["coords"])
                same = (vq > 0 and h["q"] > 0) or (vq < 0 and h["q"] < 0)
                if same and d <= args.conflict_mid:
                    sev = "STRONG" if d <= args.conflict_strong else "MID"
                    protected = vr in hs_contact_vhh
                    in_cdr = vr in cdr_res
                    buried = False
                    try:
                        buried = vhh[vr].sasa < 10.0
                    except Exception:
                        pass

                    # FIX 4: protected default = skip
                    if protected and not args.allow_protected_neutralize:
                        rows.append([name, vr, v1, h["res"], h["aa"], "same-charge", sev,
                                     round(d, 2), "", "", "", "PROTECTED(hotspot-contact): skip"])
                        continue

                    flip = ["E", "D"] if vq > 0 else ["K", "R"]
                    neut = ["Q", "N", "S", "A", "G"]
                    cands = []
                    if not (protected or in_cdr):
                        for m in flip:
                            cands.append((m, "salt-bridge", 2.0))
                    for m in neut:
                        cands.append((m, "neutralize", 1.0))

                    best = None
                    for m, why, gain in cands:
                        score = gain
                        dv = abs(VOL[m] - VOL[v1])
                        if dv > 60:
                            score -= 0.5
                        # FIX 1: hydrophobic penalty (HYDRO now used)
                        score -= max(0.0, HYDRO[m]) * 0.3
                        if buried and m in ("K", "R", "D", "E", "H"):
                            score -= 1.0
                        # FIX Bug3: protected -> heavy penalty only for flip;
                        # neutralize gets a light downgrade so the flag is usable.
                        if protected and why == "salt-bridge":
                            score -= 5.0
                        elif protected:
                            score -= 0.5
                        # 在打分中，如果冲突是 STRONG 且 <4.0 Å，G 获得"空间释放"加分
                        if sev == "STRONG" and d < 4.0 and m == "G":
                            score += 0.8

                        newq = -1 if m in ("D", "E") else (+1 if m in ("K", "R") else 0)

                        # BIOLOGY a: VHH-internal charge context
                        for g2 in gv:
                            if g2["res"] == vr:
                                continue
                            di = gdist(g["coords"], g2["coords"])
                            if di <= args.internal_cutoff:
                                qp = g2["q"]
                                if vq * qp < 0:
                                    score -= 1.0
                                if newq != 0 and newq * qp > 0:
                                    score -= 1.5
                                if newq != 0 and newq * qp < 0:
                                    score -= (1.0 if in_cdr else 0.3) #保持 VHH 构象最小扰动

                        # BIOLOGY b: new interface conflicts after mutation
                        if newq != 0:
                            for h2 in gt:
                                if h2["res"] == h["res"]:
                                    continue
                                d2 = gdist(g["coords"], h2["coords"])
                                if d2 <= args.conflict_mid and newq * h2["q"] > 0:
                                    score -= 1.5

                        if best is None or score > best[3]:
                            best = (m, why, dv, score)

                    if best is None or best[3] <= 0:
                        mut, why, sc = "", "no-safe-mutation", ""
                        note = "no-safe-mutation"
                    else:
                        mut, why, sc = best[0], best[1], round(best[3], 2)
                        note = ("PROTECTED-neutralize-only" if protected else
                                ("CDR-neutralize-only" if in_cdr else "mutate"))
                    rows.append([name, vr, v1, h["res"], h["aa"], "same-charge", sev,
                                 round(d, 2), mut, why, sc, note])
                elif (not same) and d <= args.saltbridge:
                    rows.append([name, vr, v1, h["res"], h["aa"],
                                 "salt-bridge(favorable)", "-", round(d, 2), "", "", "", "keep"])

    with open(args.out, "w") as f:
        f.write("design\tvhh_res\tvhh_aa\ttgt_res\ttgt_aa\tpair_type\tseverity\t"
                "dist_A\tsuggested_mut\tstrategy\tscore\tnote\n")
        for r in rows:
            f.write("\t".join(str(x) for x in r) + "\n")

    n_conf = sum(1 for r in rows if r[5] == "same-charge")
    n_sb = sum(1 for r in rows if "favorable" in str(r[5]))
    n_skip = sum(1 for r in rows if "skip" in str(r[11]))
    print("Conflicts: %d | Favorable salt bridges: %d | Protected-skipped: %d" % (n_conf, n_sb, n_skip))
    print("Saved: %s" % args.out)
    for r in rows:
        if r[5] != "same-charge":
            continue
        mut = f" {r[8]}" if r[8] else "  -"
        strategy = r[9] if r[9] else "-"
        score = r[10] if r[10] not in ("", None) else "-"
        print(f"  {r[0]:<14} VHH {r[2]}{r[1]} <-> TGT {r[4]}{r[3]} "
              f"{r[6]} d={r[7]:.2f} -> {mut} ({strategy}, {score}) [{r[11]}]")


if __name__ == "__main__":
    main()
