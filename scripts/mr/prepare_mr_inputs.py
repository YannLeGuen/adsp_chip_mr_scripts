#!/usr/bin/env python3
"""Prepare allele-matched, exposure-specific inputs for CHIP <-> AD MR.

CHIP input: SNP (chromosome:position:REF:ALT), Beta_unadj, SE_unadj,
P_unadj, AAF. AD input: chromosome, base_pair_location, effect_allele,
other_allele, beta, standard_error, p_value, effect_allele_frequency.

Both files must use the declared genome assembly and normalized forward-strand
alleles. The CHIP ALT allele is the effect allele.
No effects, allele labels, or frequencies are flipped here: TwoSampleMR handles
harmonization in run_mr.R. Matching requires chromosome, position AND alleles.
"""

import argparse
from pathlib import Path

import numpy as np
import pandas as pd


CHIP_COLUMNS = ["SNP", "Beta_unadj", "SE_unadj", "P_unadj", "AAF"]
AD_COLUMNS = [
    "chromosome", "base_pair_location", "effect_allele", "other_allele",
    "beta", "standard_error", "p_value", "effect_allele_frequency",
]
STANDARD_COLUMNS = ["SNP", "chr", "pos", "effect_allele", "other_allele",
                    "beta", "se", "pval", "eaf"]
DIRECTIONS = {"chip_to_ad": ("CHIP", "AD"), "ad_to_chip": ("AD", "CHIP")}


def normalize(chunk, trait):
    """Validate one input chunk; return usable rows and explicit QC counts."""
    df = pd.DataFrame(index=chunk.index)
    if trait == "CHIP":
        parts = chunk["SNP"].astype("string").str.split(":", expand=True).reindex(columns=range(4))
        df["chr"], df["pos"] = parts[0], parts[1]
        df["effect_allele"], df["other_allele"] = parts[3], parts[2]
        mapping = {"beta": "Beta_unadj", "se": "SE_unadj",
                   "pval": "P_unadj", "eaf": "AAF"}
    else:
        df["chr"], df["pos"] = chunk["chromosome"], chunk["base_pair_location"]
        df["effect_allele"] = chunk["effect_allele"]
        df["other_allele"] = chunk["other_allele"]
        mapping = {"beta": "beta", "se": "standard_error",
                   "pval": "p_value", "eaf": "effect_allele_frequency"}
    for target, source in mapping.items():
        df[target] = pd.to_numeric(chunk[source], errors="coerce")
    df["chr"] = (df["chr"].astype("string").str.strip().str.replace(r"(?i)^chr", "", regex=True)
                 .str.replace(r"\.0$", "", regex=True).str.upper().replace({"M": "MT"}))
    df["pos"] = pd.to_numeric(df["pos"], errors="coerce")
    for column in ("effect_allele", "other_allele"):
        df[column] = df[column].astype("string").str.strip().str.upper()

    valid = df["chr"].isin([str(i) for i in range(1, 23)] + ["X", "Y", "MT"])
    if trait == "CHIP":
        valid &= chunk["SNP"].astype("string").str.count(":").eq(3).fillna(False)
    valid &= np.isfinite(df["pos"]) & df["pos"].gt(0) & df["pos"].mod(1).eq(0)
    for column in ("effect_allele", "other_allele"):
        valid &= df[column].str.fullmatch("[ACGT]+", na=False)
    valid &= df["effect_allele"].ne(df["other_allele"])
    valid &= np.isfinite(df["beta"]) & np.isfinite(df["se"]) & df["se"].gt(0)
    valid &= np.isfinite(df["pval"]) & df["pval"].between(0, 1)
    df = df.loc[valid].copy()
    # Missing EAF is allowed by TwoSampleMR. Invalid supplied EAF is made missing
    # and counted, rather than deleting an otherwise usable forward-strand SNP.
    bad_eaf = df["eaf"].notna() & ~df["eaf"].between(0, 1)
    df.loc[bad_eaf, "eaf"] = np.nan
    df["pos"] = df["pos"].astype("int64")
    first = np.minimum(df["effect_allele"], df["other_allele"])
    second = np.maximum(df["effect_allele"], df["other_allele"])
    # Alleles are alphabetically sorted only in this join key, not in the data.
    # TwoSampleMR normalizes SNP identifiers to lowercase; do so here too so
    # audit joins use exactly the same identifier before and after formatting.
    df["SNP"] = (df["chr"] + ":" + df["pos"].astype(str) + ":" + first + ":" + second).str.lower()
    qc = {"rows_read": len(chunk), "valid_rows": len(df),
          "invalid_rows_removed": int((~valid).sum()),
          "invalid_eaf_set_missing": int(bad_eaf.sum()),
          "rows_p_ge_0_05": int(df["pval"].ge(0.05).sum())}
    return df[STANDARD_COLUMNS], qc


def scan_gwas(path, trait, cutoff=None, keys=None, chunksize=250_000):
    """Stream large GWAS files, keeping significant or requested rows only."""
    columns = CHIP_COLUMNS if trait == "CHIP" else AD_COLUMNS
    separator = "," if trait == "CHIP" else "\t"
    pieces, total = [], {}
    reader = pd.read_csv(path, sep=separator, usecols=columns, dtype=str,
                         chunksize=chunksize, compression="infer")
    for chunk in reader:
        df, qc = normalize(chunk, trait)
        for name, count in qc.items():
            total[name] = total.get(name, 0) + count
        keep = df["pval"].lt(cutoff) if keys is None else df["SNP"].isin(keys)
        if keep.any():
            pieces.append(df.loc[keep])
    selected = pd.concat(pieces, ignore_index=True) if pieces else pd.DataFrame(columns=STANDARD_COLUMNS)
    n_before = len(selected)
    selected = selected.drop_duplicates()
    conflicts = selected.loc[selected["SNP"].duplicated(keep=False), "SNP"].unique()
    if len(conflicts):
        example = ", ".join(conflicts[:5])
        raise ValueError(f"{trait}: conflicting summary rows for the same variant: {example}. "
                         "Resolve duplicates in the source GWAS; do not arbitrarily choose the first row.")
    total["identical_selected_duplicates_removed"] = n_before - len(selected)
    total["selected_rows"] = len(selected)
    return selected.reset_index(drop=True), total


def write_role(df, trait, role, path):
    result = df.rename(columns={c: f"{c}.{role}" for c in STANDARD_COLUMNS if c != "SNP"})
    result[role], result[f"id.{role}"] = trait, trait
    result.to_csv(path, sep="\t", index=False, na_rep="NA")


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--chip-gwas", type=Path, default=Path("/main/data/mr/chip_gwas.csv"))
    parser.add_argument("--ad-gwas", type=Path, default=Path("/main/data/mr/ad_gwas.tsv.gz"))
    parser.add_argument("--output-dir", type=Path, default=Path("/main/results/mr/inputs"))
    parser.add_argument("--direction", choices=["both", *DIRECTIONS], default="both")
    parser.add_argument("--chip-genome-wide", action="store_true",
                        help="Confirm CHIP input is genome-wide, not a table selected for CHIP significance. Required for AD -> CHIP.")
    parser.add_argument("--chip-build", choices=["GRCh37", "GRCh38"], default="GRCh38")
    parser.add_argument("--ad-build", choices=["GRCh37", "GRCh38"], default="GRCh38")
    parser.add_argument("--include-relaxed", action="store_true",
                        help="Also prepare p < 5e-6 sensitivity candidates (primary threshold is 5e-8).")
    parser.add_argument("--chunksize", type=int, default=250_000)
    args = parser.parse_args()
    directions = list(DIRECTIONS) if args.direction == "both" else [args.direction]
    if "ad_to_chip" in directions and not args.chip_genome_wide:
        parser.error("AD -> CHIP requires genome-wide CHIP summary statistics. Supply --chip-genome-wide "
                     "only after checking the input provenance, or use --direction chip_to_ad.")
    if args.chip_build != args.ad_build:
        parser.error("Genome assemblies differ. Harmonize coordinates before running; this script does not lift over.")
    if args.chunksize < 1:
        parser.error("--chunksize must be positive")
    paths = {"CHIP": args.chip_gwas, "AD": args.ad_gwas}
    for path in paths.values():
        if not path.is_file():
            parser.error(f"Required GWAS file is missing: {path}")
    cutoff = 5e-6 if args.include_relaxed else 5e-8
    candidates, qc_rows, outputs = {}, [], []
    for direction in directions:
        exposure, _ = DIRECTIONS[direction]
        print(f"Selecting {exposure} exposure candidates at p < {cutoff:g}", flush=True)
        candidates[exposure], qc = scan_gwas(paths[exposure], exposure, cutoff=cutoff,
                                             chunksize=args.chunksize)
        qc_rows.append({"trait": exposure, "pass": "exposure_selection", **qc})
        if exposure == "CHIP" and "ad_to_chip" in directions and not qc.get("rows_p_ge_0_05", 0):
            parser.error("CHIP input has no valid rows with p >= 0.05 and appears significance-filtered. "
                         "A reverse analysis requires unscreened, genome-wide CHIP outcome data.")
    # In reverse-only mode CHIP is read as an outcome in the loop below; apply
    # the same input-coverage check there, before writing any prepared files.
    for direction in directions:
        exposure, outcome = DIRECTIONS[direction]
        exp = candidates[exposure]
        print(f"Matching {outcome} outcomes to {len(exp)} {exposure} candidates", flush=True)
        out, qc = scan_gwas(paths[outcome], outcome, keys=set(exp["SNP"]), chunksize=args.chunksize)
        qc_rows.append({"trait": outcome, "pass": f"outcome_for_{exposure}", **qc})
        if outcome == "CHIP" and not qc.get("rows_p_ge_0_05", 0):
            parser.error("CHIP outcome input appears significance-filtered (no valid p >= 0.05). "
                         "Use genome-wide CHIP summary statistics for AD -> CHIP.")
        outputs.append((direction, exposure, outcome, exp, out))

    args.output_dir.mkdir(parents=True, exist_ok=True)
    manifest = []
    for direction, exposure, outcome, exp, out in outputs:
        exp_name, out_name = f"{direction}_exposure.tsv", f"{direction}_outcome.tsv"
        write_role(exp, exposure, "exposure", args.output_dir / exp_name)
        write_role(out, outcome, "outcome", args.output_dir / out_name)
        manifest.append({"direction": direction, "exposure_file": exp_name,
                         "outcome_file": out_name, "candidate_p_cutoff": cutoff,
                         "genome_build": args.chip_build, "n_exposure_candidates": len(exp),
                         "n_matched_outcomes": len(out),
                         "chip_genome_wide_confirmed": args.chip_genome_wide})
    pd.DataFrame(manifest).to_csv(args.output_dir / "analysis_inputs.tsv", sep="\t", index=False)
    pd.DataFrame(qc_rows).to_csv(args.output_dir / "input_qc.tsv", sep="\t", index=False)
    print(f"Prepared {len(manifest)} direction(s). Run run_mr.R on {args.output_dir}.")


if __name__ == "__main__":
    try:
        main()
    except (ValueError, OSError) as error:
        raise SystemExit(f"ERROR: {error}") from error
