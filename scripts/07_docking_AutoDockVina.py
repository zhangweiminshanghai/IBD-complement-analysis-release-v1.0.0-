#!/usr/bin/env python3
"""
07_docking_AutoDockVina.py

Purpose
-------
Molecular docking of nafamostat mesylate into the catalytic (S1) pocket of the
C1r/C1s serine protease domain, with complement factor I (CFI) docked as a
selectivity counter-target. Reproduces the structural-pharmacology panel of the
study:

    nafamostat - C1r/C1s catalytic domain :  dG ~ -8.8 kcal/mol
    nafamostat - CFI serine protease      :  dG ~ -6.8 kcal/mol
    selectivity window                    :  ~2 kcal/mol in favour of C1r/C1s
    catalytic triad (chymotrypsin numbering): HIS57 / ASP102 / SER195

The script (a) prepares the receptor(s), (b) prepares the ligand, (c) writes an
AutoDock Vina configuration whose search box is centred on the catalytic triad,
(d) runs Vina when the binary is available, and (e) parses the poses. When the
Vina binary (or a preparation tool) is missing the script degrades gracefully:
it still writes every input file and prints exact installation/run instructions
instead of failing.

Inputs
------
data/raw/structures/2XRC_clean.pdb        C1r/C1s serine protease domain
                                          (PDB 2XRC, cleaned: no waters/ligands)
data/raw/structures/AF-P05156-F1-model_v4.pdb
                                          Complement factor I AlphaFold model
                                          (UniProt P05156), counter-target;
                                          optional - the run skips it if absent.
                                          https://alphafold.ebi.ac.uk/entry/P05156
data/raw/ligands/nafamostat.sdf           3D conformer of nafamostat
                                          (PubChem CID 4413)

The catalytic site is located automatically: first by chymotrypsin numbering
(HIS57/ASP102/SER195), then - for structures using sequential UniProt numbering
such as 2XRC - by detecting the Ser-His-Asp charge-relay geometry, then by the
nucleophile-elbow sequence motif, and only as a last resort by blind docking.

Outputs
-------
results/07_docking/receptor_<target>.pdbqt
results/07_docking/nafamostat.pdbqt
results/07_docking/vina_config_<target>.txt
results/07_docking/vina_out_<target>.pdbqt
results/07_docking/vina_log_<target>.txt
results/07_docking/docking_scores.csv
results/07_docking/docking_summary.json
results/07_docking/binding_site_<target>.json

Figure/Table
------------
Figure 7A-C (docked pose, interaction map) ; Table 4 (binding energies and
selectivity window). Figure rendering itself is done in PyMOL/Discovery Studio
from results/07_docking/vina_out_*.pdbqt.

Key params
----------
exhaustiveness = 32, num_modes = 9, energy_range = 4 kcal/mol
box size       = 22 x 22 x 22 A centred on the catalytic triad centroid
seed           = 2024 (random_state / --seed passed to Vina)
scoring        = vina (default force field)

Approx runtime
--------------
~2-6 min per target with exhaustiveness = 32 on 8 CPU threads.
< 5 s when Vina is absent (preparation + instructions only).

Required packages / tools
-------------------------
python >= 3.9 (standard library only for the core logic)
optional python : numpy >= 1.24, pandas >= 2.0 (nicer tables),
                  meeko >= 0.5 (ligand pdbqt), rdkit >= 2023.09
external tools  : AutoDock Vina >= 1.2.5   (https://github.com/ccsb-scripps/AutoDock-Vina)
                  Open Babel >= 3.1        (obabel, receptor/ligand conversion)
                  ADFR suite               (prepare_receptor, preferred for receptors)

Author: IBD complement project
"""

from __future__ import annotations

import argparse
import csv
import json
import os
import re
import shutil
import subprocess
import sys
from dataclasses import dataclass, asdict, field
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

RANDOM_STATE = 2024

# Catalytic triad in chymotrypsinogen numbering (C1r/C1s serine protease domain)
CATALYTIC_TRIAD: Tuple[Tuple[str, int], ...] = (("HIS", 57), ("ASP", 102), ("SER", 195))

# Reference values reported in the manuscript (used for the comparison column)
REFERENCE_DG: Dict[str, float] = {"C1rC1s": -8.8, "CFI": -6.8}


# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
def find_repo_root(start: Optional[Path] = None) -> Path:
    """Walk upwards until a directory containing both data/ and scripts/ is found."""
    cur = (start or Path(__file__).resolve().parent)
    for _ in range(5):
        if (cur / "data").is_dir() and (cur / "scripts").is_dir():
            return cur
        cur = cur.parent
    return Path(__file__).resolve().parent.parent


ROOT = find_repo_root()
STRUCT_DIR = ROOT / "data" / "raw" / "structures"
LIGAND_DIR = ROOT / "data" / "raw" / "ligands"
OUT_DIR = ROOT / "results" / "07_docking"


@dataclass
class Target:
    key: str
    name: str
    receptor: Path
    chain: Optional[str] = None          # None -> auto-detect
    reference_dg: Optional[float] = None


@dataclass
class BindingSite:
    center: Tuple[float, float, float]
    size: Tuple[float, float, float]
    residues: List[str] = field(default_factory=list)
    method: str = "catalytic_triad_centroid"


# ---------------------------------------------------------------------------
# PDB parsing helpers (standard library only)
# ---------------------------------------------------------------------------
def parse_pdb_atoms(pdb_path: Path) -> List[dict]:
    """Minimal fixed-column PDB ATOM/HETATM parser."""
    atoms: List[dict] = []
    with open(pdb_path, "r", errors="ignore") as fh:
        for line in fh:
            if not line.startswith(("ATOM  ", "HETATM")):
                continue
            try:
                atoms.append(
                    {
                        "record": line[0:6].strip(),
                        "name": line[12:16].strip(),
                        "resname": line[17:20].strip().upper(),
                        "chain": line[21].strip() or "A",
                        "resseq": int(line[22:26]),
                        "x": float(line[30:38]),
                        "y": float(line[38:46]),
                        "z": float(line[46:54]),
                        "element": line[76:78].strip(),
                    }
                )
            except (ValueError, IndexError):
                continue
    if not atoms:
        raise ValueError(f"No ATOM records parsed from {pdb_path}")
    return atoms


def centroid(coords: Sequence[Tuple[float, float, float]]) -> Tuple[float, float, float]:
    n = float(len(coords))
    return (
        round(sum(c[0] for c in coords) / n, 3),
        round(sum(c[1] for c in coords) / n, 3),
        round(sum(c[2] for c in coords) / n, 3),
    )


def detect_triad_chain(atoms: List[dict]) -> Optional[str]:
    """Return the chain whose residue numbering matches the full catalytic triad."""
    chains = sorted({a["chain"] for a in atoms})
    for ch in chains:
        hits = 0
        for resname, resseq in CATALYTIC_TRIAD:
            if any(a["chain"] == ch and a["resseq"] == resseq and a["resname"] == resname
                   for a in atoms):
                hits += 1
        if hits == len(CATALYTIC_TRIAD):
            return ch
    # partial match fallback: chain with the most triad hits
    best, best_hits = None, 0
    for ch in chains:
        hits = sum(
            1
            for resname, resseq in CATALYTIC_TRIAD
            if any(a["chain"] == ch and a["resseq"] == resseq and a["resname"] == resname
                   for a in atoms)
        )
        if hits > best_hits:
            best, best_hits = ch, hits
    return best


def find_triad_by_geometry(atoms: List[dict], chain: Optional[str] = None,
                           hbond_cut: float = 3.6) -> Optional[List[dict]]:
    """
    Structure-based detection of a Ser-His-Asp(Glu) catalytic triad.

    Crystal structures such as PDB 2XRC use sequential (UniProt) numbering rather
    than chymotrypsin numbering, so HIS57/ASP102/SER195 cannot be located by
    residue number. Instead we use the defining geometry of the charge-relay
    system:  SER OG ... HIS NE2  and  HIS ND1 ... ASP OD1/OD2 (or GLU OE1/OE2),
    both within hydrogen-bonding distance.

    Returns every atom of the three triad residues, or None.
    """
    def key(a: dict) -> Tuple[str, int]:
        return (a["chain"], a["resseq"])

    def dist(a: dict, b: dict) -> float:
        return ((a["x"] - b["x"]) ** 2 + (a["y"] - b["y"]) ** 2 + (a["z"] - b["z"]) ** 2) ** 0.5

    pool = [a for a in atoms if chain is None or a["chain"] == chain]
    ser_og = [a for a in pool if a["resname"] == "SER" and a["name"] == "OG"]
    his_ne2 = [a for a in pool if a["resname"] == "HIS" and a["name"] == "NE2"]
    his_nd1 = {key(a): a for a in pool if a["resname"] == "HIS" and a["name"] == "ND1"}
    acidic = [a for a in pool if (a["resname"] == "ASP" and a["name"] in ("OD1", "OD2"))
              or (a["resname"] == "GLU" and a["name"] in ("OE1", "OE2"))]

    best: Optional[Tuple[float, Tuple[str, int], Tuple[str, int], Tuple[str, int]]] = None
    for og in ser_og:
        for ne2 in his_ne2:
            d1 = dist(og, ne2)
            if d1 > hbond_cut:
                continue
            nd1 = his_nd1.get(key(ne2))
            if nd1 is None:
                continue
            for od in acidic:
                d2 = dist(nd1, od)
                if d2 <= hbond_cut:
                    score = d1 + d2
                    if best is None or score < best[0]:
                        best = (score, key(og), key(ne2), key(od))
    if best is None:
        return None
    wanted = {best[1], best[2], best[3]}
    site = [a for a in atoms if key(a) in wanted]
    return site or None


def find_catalytic_serine_like_site(atoms: List[dict]) -> Optional[List[dict]]:
    """
    Sequence-motif fallback: locate the serine-protease nucleophile elbow
    (G-D/N-S-G-G-A/P/S/C) and return that SER plus the nearest HIS/ASP side chains.
    """
    by_chain: Dict[str, List[dict]] = {}
    for a in atoms:
        by_chain.setdefault(a["chain"], []).append(a)

    three_to_one = {
        "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q",
        "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K",
        "MET": "M", "PHE": "F", "PRO": "P", "SER": "S", "THR": "T", "TRP": "W",
        "TYR": "Y", "VAL": "V",
    }
    for chain, chain_atoms in by_chain.items():
        residues = sorted({(a["resseq"], a["resname"]) for a in chain_atoms})
        seq = "".join(three_to_one.get(r[1], "X") for r in residues)
        m = re.search(r"G[DN]SGG[APSC]", seq)
        if not m:
            continue
        ser_index = m.start() + 2
        ser_resseq = residues[ser_index][0]
        site = [a for a in chain_atoms if a["resseq"] == ser_resseq]
        ser_xyz = centroid([(a["x"], a["y"], a["z"]) for a in site])

        def dist2(a: dict) -> float:
            return (a["x"] - ser_xyz[0]) ** 2 + (a["y"] - ser_xyz[1]) ** 2 + (a["z"] - ser_xyz[2]) ** 2

        for target_res in ("HIS", "ASP"):
            near = [a for a in chain_atoms if a["resname"] == target_res and dist2(a) < 144]  # 12 A
            if near:
                nearest_resseq = min(near, key=dist2)["resseq"]
                site += [a for a in chain_atoms if a["resseq"] == nearest_resseq]
        return site
    return None


def compute_binding_site(pdb_path: Path, chain: Optional[str], box: float = 22.0) -> BindingSite:
    atoms = parse_pdb_atoms(pdb_path)
    use_chain = chain or detect_triad_chain(atoms)
    site_atoms: List[dict] = []
    if use_chain:
        for resname, resseq in CATALYTIC_TRIAD:
            site_atoms += [
                a for a in atoms
                if a["chain"] == use_chain and a["resseq"] == resseq and a["resname"] == resname
            ]
    method = "catalytic_triad_chymotrypsin_numbering"
    if len(site_atoms) < 3:
        alt = find_triad_by_geometry(atoms, use_chain)
        if alt:
            site_atoms, method = alt, "catalytic_triad_geometry_SerHisAsp"
    if len(site_atoms) < 3:
        alt = find_triad_by_geometry(atoms, None)
        if alt:
            site_atoms, method = alt, "catalytic_triad_geometry_SerHisAsp_anychain"
    if len(site_atoms) < 3:
        alt = find_catalytic_serine_like_site(atoms)
        if alt:
            site_atoms, method = alt, "nucleophile_elbow_motif"
    if len(site_atoms) < 3:
        # last resort: whole-protein centroid with a larger box (blind docking)
        prot = [(a["x"], a["y"], a["z"]) for a in atoms if a["record"] == "ATOM"]
        c = centroid(prot)
        print("  [warn] catalytic site not identified -> blind docking box (30 A)")
        return BindingSite(center=c, size=(30.0, 30.0, 30.0), residues=[], method="blind_centroid")

    c = centroid([(a["x"], a["y"], a["z"]) for a in site_atoms])
    residues = sorted({f'{a["resname"]}{a["resseq"]}:{a["chain"]}' for a in site_atoms})
    return BindingSite(center=c, size=(box, box, box), residues=residues, method=method)


# ---------------------------------------------------------------------------
# Preparation (receptor / ligand -> PDBQT)
# ---------------------------------------------------------------------------
def which(*names: str) -> Optional[str]:
    for n in names:
        p = shutil.which(n)
        if p:
            return p
    return None


def run_cmd(cmd: List[str], log_path: Optional[Path] = None, timeout: int = 3600) -> Tuple[int, str]:
    print("  $ " + " ".join(str(c) for c in cmd))
    try:
        proc = subprocess.run([str(c) for c in cmd], capture_output=True, text=True, timeout=timeout)
    except FileNotFoundError:
        return 127, f"executable not found: {cmd[0]}"
    except subprocess.TimeoutExpired:
        return 124, "timeout"
    out = (proc.stdout or "") + (proc.stderr or "")
    if log_path:
        log_path.write_text(out, encoding="utf-8")
    return proc.returncode, out


def prepare_receptor(pdb_path: Path, out_pdbqt: Path) -> bool:
    """PDB -> PDBQT via prepare_receptor (ADFR) or Open Babel."""
    if out_pdbqt.exists() and out_pdbqt.stat().st_size > 0:
        print(f"  receptor pdbqt already present: {out_pdbqt.name}")
        return True
    tool = which("prepare_receptor", "prepare_receptor4.py")
    if tool:
        rc, _ = run_cmd([tool, "-r", pdb_path, "-o", out_pdbqt, "-A", "hydrogens"])
        if rc == 0 and out_pdbqt.exists():
            return True
    obabel = which("obabel")
    if obabel:
        rc, _ = run_cmd([obabel, "-ipdb", str(pdb_path), "-opdbqt", "-O", str(out_pdbqt),
                         "-xr", "-p", "7.4", "--partialcharge", "gasteiger"])
        if rc == 0 and out_pdbqt.exists():
            return True
    print("  [warn] no receptor preparation tool found (prepare_receptor / obabel)")
    return False


def prepare_ligand(sdf_path: Path, out_pdbqt: Path) -> bool:
    """SDF -> PDBQT via Meeko (preferred) or Open Babel."""
    if out_pdbqt.exists() and out_pdbqt.stat().st_size > 0:
        print(f"  ligand pdbqt already present: {out_pdbqt.name}")
        return True
    meeko = which("mk_prepare_ligand.py")
    if meeko:
        rc, _ = run_cmd([meeko, "-i", sdf_path, "-o", out_pdbqt])
        if rc == 0 and out_pdbqt.exists():
            return True
    obabel = which("obabel")
    if obabel:
        rc, _ = run_cmd([obabel, "-isdf", str(sdf_path), "-opdbqt", "-O", str(out_pdbqt),
                         "-h", "-p", "7.4", "--partialcharge", "gasteiger"])
        if rc == 0 and out_pdbqt.exists():
            return True
    print("  [warn] no ligand preparation tool found (meeko / obabel)")
    return False


# ---------------------------------------------------------------------------
# Vina configuration and execution
# ---------------------------------------------------------------------------
def write_vina_config(cfg_path: Path, receptor: Path, ligand: Path, out_pdbqt: Path,
                      site: BindingSite, exhaustiveness: int, num_modes: int,
                      energy_range: float, cpu: int, seed: int) -> Path:
    cfg = f"""# AutoDock Vina configuration - nafamostat docking
# generated by scripts/07_docking_AutoDockVina.py (seed = {seed})
receptor = {receptor.as_posix()}
ligand   = {ligand.as_posix()}
out      = {out_pdbqt.as_posix()}

# search box centred on: {", ".join(site.residues) or site.method}
center_x = {site.center[0]}
center_y = {site.center[1]}
center_z = {site.center[2]}
size_x   = {site.size[0]}
size_y   = {site.size[1]}
size_z   = {site.size[2]}

exhaustiveness = {exhaustiveness}
num_modes      = {num_modes}
energy_range   = {energy_range}
cpu            = {cpu}
seed           = {seed}
"""
    cfg_path.write_text(cfg, encoding="utf-8")
    print(f"  wrote Vina config: {cfg_path}")
    return cfg_path


def parse_vina_output(pdbqt_path: Path) -> List[dict]:
    """Parse 'REMARK VINA RESULT: dG rmsd_lb rmsd_ub' lines from an output PDBQT."""
    poses: List[dict] = []
    if not pdbqt_path.exists():
        return poses
    with open(pdbqt_path, "r", errors="ignore") as fh:
        for line in fh:
            if line.startswith("REMARK VINA RESULT:"):
                parts = line.split()
                try:
                    poses.append(
                        {
                            "mode": len(poses) + 1,
                            "affinity_kcal_per_mol": float(parts[3]),
                            "rmsd_lb": float(parts[4]),
                            "rmsd_ub": float(parts[5]),
                        }
                    )
                except (IndexError, ValueError):
                    continue
    return poses


def parse_vina_log(log_text: str) -> List[dict]:
    """Fallback parser for the Vina stdout table."""
    poses: List[dict] = []
    for line in log_text.splitlines():
        m = re.match(r"^\s*(\d+)\s+(-?\d+\.\d+)\s+(\d+\.\d+)\s+(\d+\.\d+)\s*$", line)
        if m:
            poses.append(
                {
                    "mode": int(m.group(1)),
                    "affinity_kcal_per_mol": float(m.group(2)),
                    "rmsd_lb": float(m.group(3)),
                    "rmsd_ub": float(m.group(4)),
                }
            )
    return poses


def vina_instructions(cfg_path: Path) -> str:
    return (
        "\n".join(
            [
                "AutoDock Vina binary not found on PATH.",
                "All input files have been prepared; run the docking manually:",
                "",
                f"    vina --config {cfg_path.as_posix()}",
                "",
                "Install Vina:",
                "  * conda : conda install -c conda-forge -c bioconda autodock-vina",
                "  * binary: https://github.com/ccsb-scripps/AutoDock-Vina/releases",
                "  * macOS : brew install autodock-vina",
                "",
                "Preparation tools (if the pdbqt files are missing):",
                "  * conda install -c conda-forge openbabel   # provides obabel",
                "  * pip install meeko                        # mk_prepare_ligand.py",
                "  * ADFR suite (prepare_receptor): https://ccsb.scripps.edu/adfr/downloads/",
                "",
                "After Vina finishes, re-run this script (it parses existing outputs) or:",
                "    python scripts/docking_postprocess.py --results-dir results/07_docking",
            ]
        )
    )


def postprocess(out_dir: Path, scores_csv: Path, cutoff: float = 4.0) -> None:
    """
    Delegate pose/contact parsing to scripts/docking_postprocess.py, which takes a
    single Vina output file:
        python docking_postprocess.py <vina_out.pdbqt> --cutoff 4.0 --outdir <dir>
    and writes <stem>_docking_results.csv and <stem>_contacts.csv.
    """
    helper = Path(__file__).resolve().parent / "docking_postprocess.py"
    if not helper.exists():
        print(f"  docking_postprocess.py not present - scores written to {scores_csv.name}")
        return
    outputs = sorted(out_dir.glob("vina_out_*.pdbqt")) + sorted(out_dir.glob("vina_log_*.txt"))
    outputs = [p for p in outputs if p.stat().st_size > 0]
    if not outputs:
        print("  no Vina output files to post-process yet; run Vina first, then:")
        print(f"    python {helper.as_posix()} <vina_out.pdbqt> --outdir {out_dir.as_posix()}")
        return
    for f in outputs:
        rc, out = run_cmd([sys.executable, str(helper), str(f),
                           "--cutoff", str(cutoff), "--outdir", str(out_dir)])
        if out.strip():
            print(out.strip())
        if rc != 0:
            print(f"  [warn] docking_postprocess.py failed on {f.name} (rc = {rc})")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def dock_target(target: Target, ligand_pdbqt: Path, args) -> dict:
    print(f"\n=== {target.name} ===")
    result = {
        "target": target.key,
        "target_name": target.name,
        "receptor_pdb": str(target.receptor),
        "reference_dg": target.reference_dg,
        "status": "not_run",
        "best_affinity": None,
        "poses": [],
    }
    if not target.receptor.exists():
        print(f"  [skip] receptor not found: {target.receptor}")
        result["status"] = "missing_receptor"
        return result

    site = compute_binding_site(target.receptor, target.chain, box=args.box)
    print(f"  binding site ({site.method}): center = {site.center}, box = {site.size[0]} A")
    if site.residues:
        print("  site residues: " + ", ".join(site.residues))
    (OUT_DIR / f"binding_site_{target.key}.json").write_text(
        json.dumps(asdict(site), indent=2), encoding="utf-8")
    result["binding_site"] = asdict(site)

    receptor_pdbqt = OUT_DIR / f"receptor_{target.key}.pdbqt"
    ok_rec = prepare_receptor(target.receptor, receptor_pdbqt)
    out_pdbqt = OUT_DIR / f"vina_out_{target.key}.pdbqt"
    log_path = OUT_DIR / f"vina_log_{target.key}.txt"
    cfg_path = write_vina_config(
        OUT_DIR / f"vina_config_{target.key}.txt", receptor_pdbqt, ligand_pdbqt, out_pdbqt,
        site, args.exhaustiveness, args.num_modes, args.energy_range, args.cpu, args.seed)

    vina = which("vina", "vina.exe", "vina_1.2.5_win.exe")
    if not (vina and ok_rec and ligand_pdbqt.exists()) and not args.parse_only:
        if not vina:
            print("\n" + vina_instructions(cfg_path) + "\n")
            result["status"] = "vina_not_available"
        else:
            result["status"] = "preparation_incomplete"
    elif not args.parse_only:
        print(f"  running Vina ({vina}) ...")
        rc, out = run_cmd([vina, "--config", str(cfg_path)], log_path=log_path)
        result["status"] = "completed" if rc == 0 else f"vina_failed_rc{rc}"
        if rc != 0:
            print("  [warn] Vina failed:\n" + out[-2000:])

    poses = parse_vina_output(out_pdbqt)
    if not poses and log_path.exists():
        poses = parse_vina_log(log_path.read_text(errors="ignore"))
    if poses:
        result["poses"] = poses
        result["best_affinity"] = min(p["affinity_kcal_per_mol"] for p in poses)
        result["status"] = "parsed"
        print(f"  best affinity: {result['best_affinity']:.2f} kcal/mol "
              f"({len(poses)} poses; manuscript reference "
              f"{target.reference_dg} kcal/mol)")
    return result


def main(argv: Optional[Sequence[str]] = None) -> int:
    global OUT_DIR
    ap = argparse.ArgumentParser(description="Nafamostat docking into C1r/C1s and CFI.")
    ap.add_argument("--receptor", type=Path, default=STRUCT_DIR / "2XRC_clean.pdb",
                    help="primary receptor PDB (C1r/C1s serine protease domain)")
    ap.add_argument("--counter-receptor", type=Path,
                    default=STRUCT_DIR / "AF-P05156-F1-model_v4.pdb",
                    help="counter-target receptor PDB (CFI, UniProt P05156); skipped if absent")
    ap.add_argument("--ligand", type=Path, default=LIGAND_DIR / "nafamostat.sdf")
    ap.add_argument("--chain", default=None, help="receptor chain (default: auto-detect)")
    ap.add_argument("--box", type=float, default=22.0, help="cubic search box edge in Angstrom")
    ap.add_argument("--exhaustiveness", type=int, default=32)
    ap.add_argument("--num-modes", type=int, default=9)
    ap.add_argument("--energy-range", type=float, default=4.0)
    ap.add_argument("--cpu", type=int, default=max(1, (os.cpu_count() or 4) // 2))
    ap.add_argument("--seed", type=int, default=RANDOM_STATE)
    ap.add_argument("--out-dir", type=Path, default=OUT_DIR)
    ap.add_argument("--contact-cutoff", type=float, default=4.0,
                    help="distance cutoff (A) for ligand-residue contacts in post-processing")
    ap.add_argument("--parse-only", action="store_true",
                    help="skip execution, only parse existing Vina outputs")
    args = ap.parse_args(argv)

    OUT_DIR = args.out_dir
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    print(f"Repository root : {ROOT}")
    print(f"Output directory: {OUT_DIR}")

    # ligand ----------------------------------------------------------------
    ligand_pdbqt = OUT_DIR / "nafamostat.pdbqt"
    if args.ligand.exists():
        prepare_ligand(args.ligand, ligand_pdbqt)
    else:
        print(f"[warn] ligand not found: {args.ligand}\n"
              "       download PubChem CID 4413 (nafamostat) 3D SDF into data/raw/ligands/")

    targets = [
        Target("C1rC1s", "C1r/C1s serine protease domain (PDB 2XRC)",
               args.receptor, args.chain, REFERENCE_DG["C1rC1s"]),
        Target("CFI", "Complement factor I serine protease domain (AlphaFold)",
               args.counter_receptor, None, REFERENCE_DG["CFI"]),
    ]

    results = [dock_target(t, ligand_pdbqt, args) for t in targets]

    # scores table ----------------------------------------------------------
    scores_csv = OUT_DIR / "docking_scores.csv"
    with open(scores_csv, "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        w.writerow(["target", "target_name", "mode", "affinity_kcal_per_mol",
                    "rmsd_lb", "rmsd_ub", "reference_dg_manuscript", "status"])
        for r in results:
            if r["poses"]:
                for p in r["poses"]:
                    w.writerow([r["target"], r["target_name"], p["mode"],
                                p["affinity_kcal_per_mol"], p["rmsd_lb"], p["rmsd_ub"],
                                r["reference_dg"], r["status"]])
            else:
                w.writerow([r["target"], r["target_name"], "", "", "", "",
                            r["reference_dg"], r["status"]])

    # selectivity -----------------------------------------------------------
    best = {r["target"]: r["best_affinity"] for r in results if r["best_affinity"] is not None}
    selectivity = None
    if "C1rC1s" in best and "CFI" in best:
        selectivity = round(best["CFI"] - best["C1rC1s"], 2)
        print(f"\nSelectivity window (CFI - C1r/C1s): {selectivity:+.2f} kcal/mol "
              f"(manuscript: ~ +2.0 kcal/mol in favour of C1r/C1s)")

    summary = {
        "seed": args.seed,
        "exhaustiveness": args.exhaustiveness,
        "box_edge_angstrom": args.box,
        "catalytic_triad": [f"{n}{i}" for n, i in CATALYTIC_TRIAD],
        "ligand": str(args.ligand),
        "results": results,
        "best_affinity_kcal_per_mol": best,
        "selectivity_kcal_per_mol": selectivity,
        "reference_values": REFERENCE_DG,
        "vina_available": bool(which("vina", "vina.exe")),
    }
    (OUT_DIR / "docking_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(f"\nWrote {scores_csv.name} and docking_summary.json to {OUT_DIR}")

    postprocess(OUT_DIR, scores_csv, cutoff=args.contact_cutoff)
    return 0


if __name__ == "__main__":
    sys.exit(main())
