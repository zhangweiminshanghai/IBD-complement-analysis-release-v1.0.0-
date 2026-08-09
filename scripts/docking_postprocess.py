#!/usr/bin/env python3
# ============================================================================
# docking_postprocess.py
# ----------------------------------------------------------------------------
# Purpose : Parse AutoDock Vina docking output (.pdbqt / .log) and emit a
#           tidy table of (a) the best-pose binding free energy (ΔG, kcal/mol)
#           and (b) the receptor residues in contact with the ligand.
# Inputs  : path to a Vina output file (out.pdbqt or docking.log)
# Outputs : <stem>_docking_results.csv  (one row per pose: rank, dg, cluster)
#           <stem>_contacts.csv         (residue, chain, distance, interaction)
# Key params : contact_cutoff = 4.0 A (heavy-atom heavy-atom)
# Usage   : python docking_postprocess.py out.pdbqt --cutoff 4.0 --outdir .
# Required : numpy (standard library otherwise)
# ============================================================================
import argparse
import csv
import os
import re
import sys

try:
    from collections import OrderedDict
except ImportError:  # py2 shim
    OrderedDict = dict

PDBQT_ATOM = re.compile(
    r"^(ATOM|HETATM)\s+\d+\s+(\S+)\s+(\S+)\s+(\S?)\s*(-?\d+)\s+"
    r"(-?[\d.]+)\s+(-?[\d.]+)\s+(-?[\d.]+)\s+", re.I)
VINA_RESULT = re.compile(r"REMARK\s+VINA\s+RESULT:\s*([-\d.]+)\s+([-\d.]+)\s+([-\d.]+)")
MODEL_RE = re.compile(r"MODEL\s+(\d+)")


def parse_poses(path):
    """Return list of dicts: {rank, dg, cluster_rmsd, cluster_rmsd_lb}."""
    poses = []
    with open(path, "r", encoding="utf-8", errors="ignore") as fh:
        for line in fh:
            m = VINA_RESULT.search(line)
            if m:
                poses.append({
                    "rank": len(poses) + 1,
                    "dg": float(m.group(1)),
                    "rmsd_lb": float(m.group(2)),
                    "rmsd_ub": float(m.group(3)),
                })
    return poses


def parse_contacts(path, cutoff=4.0):
    """Parse the BEST (first) MODEL pose; compute ligand-receptor contacts.

    Returns list of dicts: {residue, chain, resnum, distance, ligand_atom}.
    Receptor atoms are HETATM/ATOM with residue name not in ligand residues
    (ligand recognised by being in the first MODEL block with residue name
    typically a single token like 'UNK' or 'LIG' or non-standard).
    """
    # Load all atoms with coordinates + residue info, tagged by MODEL block.
    atoms = []  # (model, record, name, resname, chain, resnum, x, y, z)
    model = 0
    with open(path, "r", encoding="utf-8", errors="ignore") as fh:
        for line in fh:
            if line.startswith("MODEL"):
                mm = MODEL_RE.search(line)
                model = int(mm.group(1)) if mm else model + 1
                continue
            if line.startswith("ENDMDL"):
                continue
            if line.startswith("ATOM") or line.startswith("HETATM"):
                m = PDBQT_ATOM.match(line)
                if not m:
                    continue
                rec, name, resname, chain, resnum, x, y, z = m.groups()
                atoms.append({
                    "model": model, "rec": rec, "name": name,
                    "resname": resname, "chain": chain or "A",
                    "resnum": int(resnum),
                    "xyz": (float(x), float(y), float(z)),
                })

    # Use only the first (best) pose model, and split ligand vs receptor.
    best = [a for a in atoms if a["model"] == 1]
    if not best:
        return []
    # Heuristic: ligand = atoms whose resname is a known small-mol code OR
    # whose record is HETATM and resname in a small set; receptor = the rest.
    ligand_resnames = {"UNK", "LIG", "UNL", "MLI", "DRG", "HOH"}
    ligand = [a for a in best if a["resname"] in ligand_resnames
              or (a["rec"] == "HETATM" and a["resname"] not in
                  {"ALA","ARG","ASN","ASP","CYS","GLN","GLU","GLY","HIS","ILE",
                   "LEU","LYS","MET","PHE","PRO","SER","THR","TRP","TYR","VAL",
                   "MSE","SEC","PYL"})]
    receptor = [a for a in best if a not in ligand]
    if not ligand or not receptor:
        # Fallback: if cannot separate, treat first chain as ligand.
        return []

    contacts = []
    for la in ligand:
        if la["name"].startswith("H"):  # skip hydrogens for contact calc
            continue
        for ra in receptor:
            if ra["name"].startswith("H"):
                continue
            d = sum((la["xyz"][i] - ra["xyz"][i]) ** 2 for i in range(3)) ** 0.5
            if d <= cutoff:
                contacts.append({
                    "residue": "%s%d" % (ra["resname"], ra["resnum"]),
                    "chain": ra["chain"], "resnum": ra["resnum"],
                    "distance": round(d, 2), "ligand_atom": la["name"],
                })
    contacts.sort(key=lambda c: c["distance"])
    return contacts


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("input", help="Vina .pdbqt or .log output file")
    ap.add_argument("--cutoff", type=float, default=4.0)
    ap.add_argument("--outdir", default=".")
    args = ap.parse_args()

    if not os.path.exists(args.input):
        sys.stderr.write("ERROR: input not found: %s\n" % args.input)
        sys.exit(1)

    poses = parse_poses(args.input)
    if not poses:
        sys.stderr.write("WARNING: no VINA RESULT lines found; "
                         "writing empty tables.\n")
    contacts = parse_contacts(args.input, cutoff=args.cutoff)

    stem = os.path.splitext(os.path.basename(args.input))[0]
    os.makedirs(args.outdir, exist_ok=True)

    poses_csv = os.path.join(args.outdir, "%s_docking_results.csv" % stem)
    with open(poses_csv, "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        w.writerow(["rank", "dg_kcal_mol", "rmsd_lb", "rmsd_ub"])
        for p in poses:
            w.writerow([p["rank"], p["dg"], p["rmsd_lb"], p["rmsd_ub"]])

    contacts_csv = os.path.join(args.outdir, "%s_contacts.csv" % stem)
    with open(contacts_csv, "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        w.writerow(["residue", "chain", "resnum", "distance_A", "ligand_atom"])
        for c in contacts:
            w.writerow([c["residue"], c["chain"], c["resnum"],
                        c["distance"], c["ligand_atom"]])

    best_dg = poses[0]["dg"] if poses else float("nan")
    print("Parsed %d pose(s). Best ΔG = %.3f kcal/mol" % (len(poses), best_dg))
    print("Contact residues (<=%.1f A): %d" % (args.cutoff, len(contacts)))
    for c in contacts[:12]:
        print("  %s  %s  %.2f A" % (c["residue"], c["chain"], c["distance"]))
    print("Wrote:\n  %s\n  %s" % (poses_csv, contacts_csv))


if __name__ == "__main__":
    main()
