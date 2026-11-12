#!/usr/bin/env bash
# Example: screen NK2R-targeting VHH designs (15 NKA-activation hotspots),
# then optimize interface charge conflicts on the survivors.
#
# Adjust --pdb_dir / chain ids / hotspots to your own target.
set -euo pipefail

PDB_DIR="${1:-predictions/}"          # dir of 2-chain PDBs (A=VHH, B=NK2R)

# NK2R activation-pocket hotspots (UniProt P21452 numbering), side-chain atoms
HOTSPOTS="86:ND2,90:ND2,93:OH,94:CB,96:NE2,114:CG1/CG2/CD1,117:SD/CE,175:OD1/OD2,180:NZ,181:SG,269:OH,270:CZ/CE1/CE2,280:OH,289:OH,293:CZ/CE1/CE2"

echo "== Step 1: hotspot coverage + terminal-burial screen =="
python scripts/hotspot_terminal_filter.py \
    --pdb_dir  "$PDB_DIR" \
    --chain_vhh A --chain_tgt B \
    --hotspots "$HOTSPOTS" \
    --cutoff 4.5 \
    --term_len 25 --term_max_hs 2 \
    --out ranking.tsv

echo
echo "== Step 2: interface charge-conflict optimization =="
python scripts/charge_conflict_optimize.py \
    --pdb_dir  "$PDB_DIR" \
    --chain_vhh A --chain_tgt B \
    --hotspots 86,90,93,94,96,114,117,175,180,181,269,270,280,289,293 \
    --iface_cutoff 5.0 \
    --allow_protected_neutralize \
    --out charge_optimization.tsv

echo
echo "Outputs: ranking.tsv, charge_optimization.tsv"
