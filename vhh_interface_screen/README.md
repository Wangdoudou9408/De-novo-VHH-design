# vhh-interface-screen

Tools for screening and locally optimizing **de-novo designed binders**
(e.g. VHH / nanobodies) against a structured target, using predicted
complex structures (AlphaFold / Boltz / RFdiffusion outputs).

Two complementary stages:

1. **`hotspot_terminal_filter.py`** — rank complexes by how many predefined
   target *hotspot* residues the binder contacts, and **reject** designs whose
   binder **N- or C-terminus buries into the hotspot pocket** (affinity tags or
   Fc fusions attached at the termini would sterically block binding).
2. **`charge_conflict_optimize.py`** — for designs that pass, detect
   **same-charge electrostatic repulsions** at the interface and propose
   single-point mutations on the **binder side only** to relieve them
   (neutralize, or flip to form a salt bridge), while protecting
   hotspot-contacting and CDR residues.

Neither tool is a black-box scorer: both emit per-residue, interpretable
tables you can inspect and act on.

---

## Install

```bash
git clone <this-repo> && cd vhh-interface-screen
python -m venv .venv && source .venv/bin/activate
pip install -r requirements.txt
# optional, for CDR/FR-aware charge optimization:
pip install anarci
```

Requires Python ≥ 3.9.

---

## Input format

A directory of **two-chain PDB files**:

| chain | content |
|-------|---------|
| `--chain_vhh` (default `A`) | the designed binder (VHH / nanobody / miniprotein) |
| `--chain_tgt` (default `B`) | the target (receptor / antigen) |

Residue numbering must be consistent with the hotspot list you supply
(see below). Structures typically come from AlphaFold-Multimer, Boltz, or
RFdiffusion+ProteinMPNN pipelines.

---

## Usage

### 1. Hotspot + terminal-burial screen

```bash
python scripts/hotspot_terminal_filter.py \
    --pdb_dir  predictions/ \
    --chain_vhh A --chain_tgt B \
    --hotspots "86:ND2,90:ND2,93:OH,94:CB,96:NE2,114:CG1/CG2/CD1,117:SD/CE,175:OD1/OD2,180:NZ,181:SG,269:OH,270:CZ/CE1/CE2,280:OH,289:OH,293:CZ/CE1/CE2" \
    --cutoff 4.5 \
    --term_len 25 --term_max_hs 2 \
    --out ranking.tsv
```

Hotspot syntax: `resnum` (all heavy atoms) or `resnum:atom1/atom2` (specific
side-chain atoms). The example above is a 15-residue GPCR activation-pocket
set; supply your own for other targets.

Output columns:
`rank, name, vhh_len, n_hotspots, hotspot_list, n_contacts, nterm_hs, cterm_hs, terminal_pass, fail_reason`

Survivors (`terminal_pass=PASS`) are ranked first by `n_hotspots`, then
`n_contacts`.

### 2. Interface charge-conflict optimization

```bash
python scripts/charge_conflict_optimize.py \
    --pdb_dir  predictions/ \
    --chain_vhh A --chain_tgt B \
    --hotspots 86,90,93,94,96,114,117,175,180,181,269,270,280,289,293 \
    --iface_cutoff 5.0 \
    --out charge_optimization.tsv
```

Useful flags:

| flag | effect |
|------|--------|
| `--allow_protected_neutralize` | hotspot-contacting residues may be *neutralized* (never flipped); default = skip them |
| `--use_anarci` | number the binder with ANARCI (IMGT); CDR positions only get neutralize |
| `--internal_cutoff` | radius for binder-internal charge-context checks (default 6 Å) |

Output columns:
`design, vhh_res, vhh_aa, tgt_res, tgt_aa, pair_type, severity, dist_A, suggested_mut, strategy, score, note`

`note` values: `mutate`, `PROTECTED-neutralize-only`, `CDR-neutralize-only`,
`PROTECTED(hotspot-contact): skip`, `no-safe-mutation`, `keep` (favourable
salt bridge).

---

## Recommended workflow

```
predicted complexes (PDB)
        │
        ▼
hotspot_terminal_filter.py      →  rank by hotspot coverage,
        │                            drop terminus-buried designs
        ▼
charge_conflict_optimize.py     →  list interface charge repulsions +
        │                            suggested binder mutations
        ▼
apply mutations in silico  →  re-predict (AF/Boltz)  →  validate (MD / assay)
```

Treat the suggested mutations as **hypotheses**: always re-predict the mutant
complex and re-check hotspot coverage before committing to wet-lab work.

---

## Notes & caveats

- Distances are computed on **heavy atoms** with a KD-tree; `query_ball_point`
  is used so *all* contacting residues are captured (not just the nearest).
- His is treated as weakly positive (+0.5) since its protonation is
  pH-dependent; adjust in `charged_groups()` if you model a specific pH.
- The terminal-burial rule is a heuristic against tag/Fc interference; if your
  construct places tags internally or uses a tag-free expression, relax
  `--term_max_hs`.
- This code does **not** run structure prediction itself; feed it predictions
  from AlphaFold / Boltz / RoseTTAFold / RFdiffusion pipelines.

## License

MIT — see [LICENSE](LICENSE).
