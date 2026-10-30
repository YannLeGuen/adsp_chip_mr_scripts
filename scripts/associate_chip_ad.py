#!/usr/bin/env python3
"""Fit the essential ADSP CHIP–AD models and write tidy regression tables.

WGS: AD ~ CHIP + sex + age + APOE2 + APOE4 + PC1 + PC2 + PC3 + PC4.
WES: European participants, age, sex, APOE dosages, and two European PCs.
Within APOE strata, APOE dosages are omitted because they are constant.
Input phenotypes must be harmonized to the documented column definitions.
"""

from __future__ import annotations

import argparse
import json
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd


def read_table(path: Path) -> pd.DataFrame:
    if not path.is_file():
        raise FileNotFoundError(f"Required input is missing: {path}")
    return pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)


def prepare_cohort(phenotypes: pd.DataFrame, status: pd.DataFrame,
                   cohort: str) -> tuple[pd.DataFrame, dict]:
    """Validate IDs/coding, join measured CHIP status, and count exclusions."""
    pcs = [f"PC{i}" for i in range(1, 5 if cohort == "WGS" else 3)]
    required = {"sample_id", "participant_id", "AD", "sex", "age", "APOE", *pcs}
    if cohort == "WES":
        required.add("ancestry")
    absent = sorted(required - set(phenotypes))
    if absent:
        raise ValueError(f"Phenotypes lack columns: {', '.join(absent)}")
    if not {"sample_id", "participant_id", "CHIP"} <= set(status):
        raise ValueError("Status file must contain sample_id, participant_id, and CHIP.")
    for label, table in [("phenotypes", phenotypes), ("status", status)]:
        if not table.sample_id.eq(table.sample_id.str.strip()).all():
            raise ValueError(f"{label} sample_id values contain surrounding whitespace; normalize IDs first.")
        if table.sample_id.str.strip().eq("").any() or table.sample_id.duplicated().any():
            raise ValueError(f"{label} must have a nonempty, unique sample_id per row.")
    if not status.participant_id.eq(status.participant_id.str.strip()).all():
        raise ValueError("Status participant_id values contain surrounding whitespace.")
    if status.participant_id.eq("").any() or status.participant_id.duplicated().any():
        raise ValueError("The status file must contain one sample per unique participant.")
    # Site QC includes all eligible blood participants. Association uses the
    # clinically annotated subset, keeping both exclusions explicit in counts.
    frame = phenotypes[sorted(required)].merge(
        status[["sample_id", "participant_id", "CHIP"]], on="sample_id",
        suffixes=("", "_status"), validate="one_to_one")
    if not frame.participant_id.eq(frame.participant_id.str.strip()).all():
        raise ValueError("participant_id values contain surrounding whitespace; normalize IDs first.")
    if frame.participant_id.str.strip().eq("").any() or frame.participant_id.duplicated().any():
        raise ValueError("Select one blood sample per participant before analysis; participant_id must be unique.")
    if not frame.participant_id.eq(frame.participant_id_status).all():
        raise ValueError("Phenotype participant_id does not match the blood-sample manifest for one or more samples.")
    frame = frame.drop(columns="participant_id_status")
    report = {"phenotype_rows": len(phenotypes), "status_rows": len(status),
              "matched_samples": len(frame),
              "status_samples_without_phenotype": len(set(status.sample_id) - set(phenotypes.sample_id)),
              "phenotype_samples_outside_status": len(set(phenotypes.sample_id) - set(status.sample_id))}
    if cohort == "WES":
        frame = frame.loc[frame.ancestry.str.strip().str.upper().eq("EUR")].copy()
    report["ancestry_eligible"] = len(frame)
    missing_tokens = {"", "NA", "NaN", "nan", "NAN", ".", "-9"}
    for column in ["AD", "sex", "age", "CHIP", *pcs]:
        source = frame[column].astype(str).str.strip()
        if column == "age":
            source = source.replace("90+", "90")
        source = source.mask(source.isin(missing_tokens))
        numeric = pd.to_numeric(source, errors="coerce")
        if (source.notna() & numeric.isna()).any():
            raise ValueError(f"{column} contains unrecognized nonnumeric values.")
        if np.isinf(numeric).any():
            raise ValueError(f"{column} contains infinite values.")
        frame[column] = numeric
    for column in ["AD", "sex", "CHIP"]:
        if not frame[column].dropna().isin([0, 1]).all():
            raise ValueError(f"{column} must be coded 0/1 (or missing); see README.")
    if (frame.age.dropna() <= 0).any():
        raise ValueError("Age must be positive or missing.")
    apoe = frame.APOE.astype(str).str.strip().str.upper().str.removeprefix("E")
    apoe = apoe.str.replace(r"\.0$", "", regex=True).mask(apoe.isin(missing_tokens))
    if not apoe.dropna().isin(["22", "23", "24", "33", "34", "44"]).all():
        raise ValueError("APOE must be 22, 23, 24, 33, 34, or 44 (or missing).")
    frame["APOE"] = apoe
    frame["APOE2"] = apoe.str.count("2").astype(float)
    frame["APOE4"] = apoe.str.count("4").astype(float)
    # The pooled and stratified models use the same complete-case cohort.
    complete = ["AD", "sex", "age", "CHIP", "APOE", *pcs]
    report["missing_CHIP_status"] = int(frame.CHIP.isna().sum())
    report["incomplete_covariates_or_status"] = int(frame[complete].isna().any(axis=1).sum())
    frame = frame.dropna(subset=complete).copy()
    report["complete_cases"] = len(frame)
    if frame.empty:
        raise ValueError("No complete, eligible participants remain.")
    return frame, report


def fit_model(frame: pd.DataFrame, label: str, covariates: list[str]) -> dict:
    """Fit an ordinary binomial GLM; never substitute a penalized estimator."""
    import statsmodels.api as sm
    import statsmodels.formula.api as smf
    from statsmodels.tools.sm_exceptions import PerfectSeparationError

    formula = "AD ~ CHIP + " + " + ".join(covariates)
    row = {"APOE": label, "N": len(frame), "N_AD": int(frame.AD.sum()),
           "N_control": int((frame.AD == 0).sum()), "N_CHIP": int(frame.CHIP.sum()),
           "formula": formula, "status": "not_fitted"}
    cells = pd.crosstab(frame.AD, frame.CHIP).reindex(index=[0, 1], columns=[0, 1], fill_value=0)
    if (cells == 0).any().any():
        row["status"] = "not_estimable_empty_case_control_CHIP_cell"
        return row
    try:
        model = smf.glm(formula, data=frame, family=sm.families.Binomial(), missing="raise")
        if np.linalg.matrix_rank(model.exog) < model.exog.shape[1]:
            row["status"] = "not_estimable_rank_deficient_design"
            return row
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            result = model.fit(maxiter=200)
        fitting_warnings = [str(item.message) for item in caught]
        if not result.converged or any("separation" in message.lower() for message in fitting_warnings):
            row["status"] = "not_estimable_nonconvergence_or_separation"
            return row
        beta, se, p = [float(value) for value in (result.params["CHIP"], result.bse["CHIP"], result.pvalues["CHIP"])]
        low, high = map(float, result.conf_int().loc["CHIP"])
        if not np.isfinite([beta, se, p, low, high]).all() or se <= 0:
            row["status"] = "not_estimable_nonfinite_estimate"
            return row
        with np.errstate(over="raise"):
            odds, ci_low, ci_high = np.exp([beta, low, high])
        row.update(beta=beta, se=se, OR=float(odds), CI_lower=float(ci_low),
                   CI_upper=float(ci_high), p_value=p, status="ok",
                   warnings="; ".join(fitting_warnings))
    except (ValueError, np.linalg.LinAlgError, FloatingPointError, PerfectSeparationError) as error:
        row["status"] = f"not_estimable_{type(error).__name__}"
    return row


def demographics(frame: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for (ad, chip), group in frame.groupby(["AD", "CHIP"], sort=True):
        rows.append({"AD": int(ad), "CHIP": int(chip), "N": len(group),
                     "age_mean": group.age.mean(), "age_SD": group.age.std(ddof=1),
                     "age_median": group.age.median(), "age_min": group.age.min(),
                     "age_max": group.age.max(), "female_N": int(group.sex.sum()),
                     **{f"APOE_{genotype}_N": int(group.APOE.eq(genotype).sum())
                        for genotype in ["22", "23", "24", "33", "34", "44"]}})
    return pd.DataFrame(rows)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cohort", choices=["WGS", "WES"], required=True)
    parser.add_argument("--phenotypes", type=Path, help="Default /main/data/adsp_{cohort}/phenotypes.tsv")
    parser.add_argument("--status", type=Path, help="Default /main/results/adsp_{cohort}/status.tsv")
    parser.add_argument("--output-dir", type=Path, help="Default /main/results/adsp_{cohort}/association")
    args = parser.parse_args()
    cohort = args.cohort.lower()
    phenotypes = args.phenotypes or Path(f"/main/data/adsp_{cohort}/phenotypes.tsv")
    status = args.status or Path(f"/main/results/adsp_{cohort}/status.tsv")
    output = args.output_dir or Path(f"/main/results/adsp_{cohort}/association")
    try:
        frame, report = prepare_cohort(read_table(phenotypes), read_table(status), args.cohort)
        pcs = [f"PC{i}" for i in range(1, 5 if args.cohort == "WGS" else 3)]
        covariates = ["sex", "age", *pcs]
        rows = [fit_model(frame, "All", [*covariates, "APOE2", "APOE4"])]
        for genotype in ["23", "33", "34", "44"]:
            rows.append(fit_model(frame.loc[frame.APOE.eq(genotype)], f"E{genotype}", covariates))
        import statsmodels
        report.update(cohort=args.cohort, age_adjusted=True,
                      software={"python": sys.version.split()[0], "pandas": pd.__version__,
                                "numpy": np.__version__, "statsmodels": statsmodels.__version__})
        output.mkdir(parents=True, exist_ok=True)
        result_columns = ["APOE", "N", "N_AD", "N_control", "N_CHIP", "beta", "se", "OR",
                          "CI_lower", "CI_upper", "p_value", "formula", "status", "warnings"]
        pd.DataFrame(rows).reindex(columns=result_columns).to_csv(
            output / "chip_ad_association.tsv", sep="\t", index=False, na_rep="NA")
        demographics(frame).to_csv(output / "demographics.tsv", sep="\t", index=False, na_rep="NA")
        (output / "analysis_counts.json").write_text(json.dumps(report, indent=2) + "\n")
        print(f"Analyzed {len(frame)} participants; results: {output}")
    except (ValueError, OSError, ImportError) as error:
        parser.exit(2, f"ERROR: {error}\n")


if __name__ == "__main__":
    main()
