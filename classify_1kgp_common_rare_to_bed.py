#!/usr/bin/env python3

import pandas as pd
from pathlib import Path

POPS = ["AFR", "AMR", "EAS", "EUR", "SAS"]

AF_GLOBAL_DIR = Path("AF_global")
AF_SUPERPOP_DIR = Path("AF_superpop")
BIM_DIR = Path("ld_pruned")

OUT_COMMON = Path("1KGP_common_LDpruned_WG.bed")
OUT_RARE = Path("1KGP_rare_LDpruned_WG.bed")
OUT_INTERMEDIATE = Path("1KGP_intermediate_LDpruned_WG.bed")
OUT_SUMMARY = Path("1KGP_common_rare_summary.tsv")

common_rows = []
rare_rows = []
intermediate_rows = []
summary_rows = []

def chrom_to_bed_chrom(x):
    x = str(x)
    return x if x.startswith("chr") else "chr" + x

for chr_num in range(1, 23):
    print(f"Processing chr{chr_num}")

    bim_file = BIM_DIR / f"chr{chr_num}.LDpruned.bim"
    global_af_file = AF_GLOBAL_DIR / f"chr{chr_num}.afreq"

    bim = pd.read_csv(
        bim_file,
        sep=r"\s+",
        header=None,
        names=["CHR", "ID", "CM", "POS", "A1", "A2"]
    )

    df = bim[["CHR", "ID", "POS"]].copy()

    # Global AF
    g_af = pd.read_csv(global_af_file, sep=r"\s+")
    g_af = g_af[["ID", "ALT_FREQS", "OBS_CT"]].rename(
        columns={"ALT_FREQS": "AF_GLOBAL", "OBS_CT": "OBS_CT_GLOBAL"}
    )
    df = df.merge(g_af, on="ID", how="left")

    # Superpopulation AF
    for pop in POPS:
        af_file = AF_SUPERPOP_DIR / f"chr{chr_num}.{pop}.afreq"

        af = pd.read_csv(af_file, sep=r"\s+")
        af = af[["ID", "ALT_FREQS", "OBS_CT"]].rename(
            columns={
                "ALT_FREQS": f"AF_{pop}",
                "OBS_CT": f"OBS_CT_{pop}"
            }
        )

        df = df.merge(af, on="ID", how="left")

    af_cols = [f"AF_{p}" for p in POPS]

    # Safety check
    if df["AF_GLOBAL"].isna().any():
        n_missing = df["AF_GLOBAL"].isna().sum()
        raise ValueError(f"Missing global AF for {n_missing} variants in chr{chr_num}")

    df[af_cols] = df[af_cols].fillna(0)

    df["n_superpop_AF_ge_1pct"] = (df[af_cols] >= 0.01).sum(axis=1)

    # Final definitions, matching thesis logic
    df["is_common"] = (df["AF_GLOBAL"] >= 0.01) & (df["n_superpop_AF_ge_1pct"] >= 2)
    df["is_rare"] = df["AF_GLOBAL"] < 0.01
    df["is_intermediate"] = ~(df["is_common"] | df["is_rare"])

    common = df[df["is_common"]].copy()
    rare = df[df["is_rare"]].copy()
    intermediate = df[df["is_intermediate"]].copy()

    def make_bed(subset):
        return pd.DataFrame({
            "chrom": subset["CHR"].apply(chrom_to_bed_chrom),
            "start": subset["POS"].astype(int) - 1,
            "end": subset["POS"].astype(int),
            "id": subset["ID"],
            "AF_GLOBAL": subset["AF_GLOBAL"],
            "n_superpop_AF_ge_1pct": subset["n_superpop_AF_ge_1pct"]
        })

    common_rows.append(make_bed(common))
    rare_rows.append(make_bed(rare))
    intermediate_rows.append(make_bed(intermediate))

    summary_rows.append({
        "chromosome": chr_num,
        "LD_pruned_total": len(df),
        "common_AFglobal_ge1pct_and_ge2superpops": len(common),
        "rare_AFglobal_lt1pct": len(rare),
        "intermediate_not_common_not_rare": len(intermediate)
    })

common_all = pd.concat(common_rows, ignore_index=True)
rare_all = pd.concat(rare_rows, ignore_index=True)
intermediate_all = pd.concat(intermediate_rows, ignore_index=True)
summary = pd.DataFrame(summary_rows)

common_all = common_all.sort_values(["chrom", "start", "end"])
rare_all = rare_all.sort_values(["chrom", "start", "end"])
intermediate_all = intermediate_all.sort_values(["chrom", "start", "end"])

common_all.to_csv(OUT_COMMON, sep="\t", header=False, index=False)
rare_all.to_csv(OUT_RARE, sep="\t", header=False, index=False)
intermediate_all.to_csv(OUT_INTERMEDIATE, sep="\t", header=False, index=False)
summary.to_csv(OUT_SUMMARY, sep="\t", index=False)

print("\nDONE")
print(f"Common BED:       {OUT_COMMON} ({len(common_all):,} variants)")
print(f"Rare BED:         {OUT_RARE} ({len(rare_all):,} variants)")
print(f"Intermediate BED: {OUT_INTERMEDIATE} ({len(intermediate_all):,} variants)")
print(f"Summary:          {OUT_SUMMARY}")

print("\nTotal check:")
print(f"Common + Rare + Intermediate = {len(common_all) + len(rare_all) + len(intermediate_all):,}")
print(f"Total LD-pruned variants     = {summary['LD_pruned_total'].sum():,}")
