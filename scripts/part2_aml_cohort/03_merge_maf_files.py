#!/usr/bin/env python3
"""
03_merge_maf_files.py
Merge all per-sample TCGA-LAML MAF files into a single MAF for maftools.

Usage: python3 scripts/part2_aml_cohort/03_merge_maf_files.py
"""

import gzip
import glob
import os

INDIR   = "data/tcga_maf"
OUTDIR  = "data/processed"
OUTFILE = os.path.join(OUTDIR, "TCGA-LAML-merged.maf")

os.makedirs(OUTDIR, exist_ok=True)

maf_files = sorted(glob.glob(os.path.join(INDIR, "*", "*.maf.gz")))
print(f"[INFO] Found {len(maf_files)} MAF files")

if not maf_files:
    raise SystemExit("[ERROR] No MAF files found")

header_written = False
total_records  = 0

with open(OUTFILE, "w") as out:
    for i, path in enumerate(maf_files):
        with gzip.open(path, "rt") as f:
            lines = f.readlines()

        comment_lines = [l for l in lines if l.startswith("#")]
        header_line   = next((l for l in lines if l.startswith("Hugo_Symbol")), None)
        data_lines    = [l for l in lines
                         if not l.startswith("#") and not l.startswith("Hugo_Symbol")]

        if not header_written:
            for c in comment_lines:
                out.write(c)
            if header_line:
                out.write(header_line)
            header_written = True

        for line in data_lines:
            out.write(line)
        total_records += len(data_lines)

        if (i + 1) % 20 == 0:
            print(f"[INFO] Processed {i+1}/{len(maf_files)} files...")

print(f"[INFO] Done. Merged MAF: {OUTFILE}")
print(f"[INFO] Total mutation records: {total_records}")