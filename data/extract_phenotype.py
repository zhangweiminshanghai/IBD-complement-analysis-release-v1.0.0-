#!/usr/bin/env python3
"""Extract CD/Control phenotype labels from a GEO series matrix (.txt).

GEO series matrices are TRANSPOSED: each !Sample_* metadata line holds
tab-separated values, one per sample (column). The disease group is taken
from !Sample_source_name_ch1 (clearest signal: 'control' / 'CDc' / 'CDi' /
'UC'). Output: data/processed/<GSE>_phenotype.csv (sample, group).
Group logic (mirrors GEO_CD.R): Crohn's (CDc/CDi) -> CD; ulcerative colitis
-> UC; control/healthy -> Control; else Other (excluded from CD-vs-Control).
"""
import csv, os, re

def grab(line):
    # drop the !Sample_xxx tag, strip surrounding quotes
    return [x.strip().strip('"') for x in line.rstrip("\n").split("\t")[1:]]

def classify(text):
    t = text.lower()
    if re.search(r"\bcd[ci]\b", t) or re.search(r"\bcd\b", t) or "crohn" in t:
        return "CD"
    if "ulcerative colitis" in t or re.search(r"\buc\b", t):
        return "UC"
    if "control" in t or "healthy" in t or "normal" in t:
        return "Control"
    return "Other"

def extract(series_matrix, out_csv):
    acc = src = None
    with open(series_matrix, "r", encoding="utf-8", errors="ignore") as fh:
        for line in fh:
            if line.startswith("!Sample_geo_accession"):
                acc = grab(line)
            elif line.startswith("!Sample_source_name_ch1"):
                src = grab(line)
    if acc is None or src is None:
        raise SystemExit("Could not find sample/metadata lines in " + series_matrix)
    n = min(len(acc), len(src))
    out = []
    for i in range(n):
        out.append((acc[i], classify(src[i])))
    with open(out_csv, "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        w.writerow(["sample", "group"])
        for s, g in out:
            w.writerow([s, g])
    print("Wrote %s (%d samples)" % (out_csv, len(out)))

if __name__ == "__main__":
    base = r"D:\IBD_Project"
    out = r"D:\IBD_Project\IBD-complement-analysis\data\processed"
    os.makedirs(out, exist_ok=True)
    extract(os.path.join(base, "GSE16879_series_matrix.txt"),
            os.path.join(out, "GSE16879_phenotype.csv"))
    extract(os.path.join(base, "GSE75214_series_matrix.txt"),
            os.path.join(out, "GSE75214_phenotype.csv"))
