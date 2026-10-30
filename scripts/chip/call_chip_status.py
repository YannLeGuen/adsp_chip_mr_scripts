#!/usr/bin/env python3
"""Extract allele depths and determine CHIP status in ADSP WGS or WES.

Variant QC is evaluated in unique blood-derived participants before clinical
or covariate exclusions. The same participants define the WES prevalence
denominator. Individual calls require VAF >2% and total allele depth >15 for
WGS or >30 for WES. Missing allele depths remain distinct from negative calls.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import math
import shutil
import statistics
import subprocess
import tempfile
from contextlib import contextmanager
from pathlib import Path

from select_chip_variants import canonical_variant, read_table


def load_samples(path: Path) -> dict[str, str]:
    """Read and validate the full eligible blood cohort before phenotype filtering."""
    frame = read_table(path)
    required = {"sample_id", "participant_id", "body_site", "technical_replicate"}
    absent = sorted(required - set(frame.columns))
    if absent:
        raise ValueError(f"Sample table lacks columns: {', '.join(absent)}")
    for column in ("sample_id", "participant_id"):
        if not frame[column].eq(frame[column].str.strip()).all():
            raise ValueError(f"{column} contains surrounding whitespace; normalize identifiers first.")
    samples = frame["sample_id"].str.strip().tolist()
    if not samples or any(not value for value in samples) or len(set(samples)) != len(samples):
        raise ValueError("sample_id values must be nonempty and unique.")
    participants = frame["participant_id"]
    if (participants == "").any() or participants.duplicated().any():
        raise ValueError("Select exactly one sample per participant before running CHIP analysis.")
    if not frame["body_site"].str.strip().str.lower().eq("blood").all():
        raise ValueError("The sample table must contain blood-derived samples only.")
    if not frame["technical_replicate"].str.strip().str.lower().isin(["0", "0.0", "false"]).all():
        raise ValueError("Exclude technical replicates from the sample table.")
    return dict(zip(samples, participants))


def load_candidates(path: Path) -> dict[tuple, dict]:
    table = read_table(path)
    required = {"SNP", "SYMBOL", "Consequence", "HGVSp", "driver_rule"}
    if required - set(table.columns):
        raise ValueError("Candidate table must be the output of select_chip_variants.py.")
    candidates = {}
    for row in table.to_dict("records"):
        key = canonical_variant(row["SNP"])
        if key in candidates:
            raise ValueError(f"Duplicate candidate variant: {row['SNP']}")
        candidates[key] = row
    if not candidates:
        raise ValueError("The candidate table is empty.")
    return candidates


def open_text(path: Path):
    return gzip.open(path, "rt") if path.suffix in {".gz", ".bgz"} else path.open()


def parse_header(lines) -> tuple[list[str], dict[str, str]]:
    contigs = {}
    for line in lines:
        if line.startswith("##contig=<ID="):
            raw = line.split("##contig=<ID=", 1)[1].split(",", 1)[0].split(">", 1)[0]
            canonical = canonical_variant(f"{raw}:1:A:C")[0]
            if canonical in contigs and contigs[canonical] != raw:
                raise ValueError(f"Ambiguous chromosome aliases in VCF header: {raw}")
            contigs[canonical] = raw
        if line.startswith("#CHROM\t"):
            fields = line.rstrip("\r\n").split("\t")
            if fields[:9] != ["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT"]:
                raise ValueError("VCF must contain FORMAT and sample columns.")
            samples = fields[9:]
            if not samples or any(not sample for sample in samples) or len(samples) != len(set(samples)):
                raise ValueError("VCF sample names must be nonempty and unique.")
            return samples, contigs
    raise ValueError("VCF #CHROM header was not found.")


@contextmanager
def vcf_stream(path: Path, candidates: dict, reader: str, scratch: Path):
    """Use indexed extraction for compressed VCF/BCF, or explicit sequential scan."""
    compressed = path.suffix in {".gz", ".bgz", ".bcf"}
    use_bcftools = reader == "bcftools" or (reader == "auto" and compressed)
    if not use_bcftools:
        if path.suffix == ".bcf":
            raise ValueError("BCF requires --reader bcftools.")
        with open_text(path) as handle:
            samples, _ = parse_header(handle)
            yield samples, handle
        return
    executable = shutil.which("bcftools")
    if executable is None:
        raise ValueError("Install bcftools for indexed VCF extraction, or explicitly use --reader scan for VCF text.")
    header = subprocess.run([executable, "view", "-h", str(path)], text=True, capture_output=True, check=True)
    samples, contigs = parse_header(header.stdout.splitlines())
    if not contigs:
        raise ValueError(f"Indexed extraction requires ##contig declarations: {path.name}")
    region_file = scratch / "query_regions.tsv"
    regions = {(contigs[key[0]], key[1]) for key in candidates if key[0] in contigs}
    if not regions:
        yield samples, iter(())
        return
    with region_file.open("w") as output:
        for chrom, pos in sorted(regions):
            output.write(f"{chrom}\t{pos}\n")
    # A .tsv regions file is 1-based. --regions-overlap 0 selects record positions.
    command = [executable, "view", "--regions-file", str(region_file), "--regions-overlap", "0", "-Ov", str(path)]
    with tempfile.TemporaryFile(mode="w+t", dir=scratch) as errors:
        process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=errors, text=True)
        try:
            assert process.stdout is not None
            yield samples, process.stdout
            process.stdout.close()
            returncode = process.wait()
            if returncode:
                errors.seek(0)
                raise ValueError(f"bcftools extraction failed for {path.name}: {errors.read().strip()}")
        finally:
            if process.poll() is None:
                process.terminate()
                process.wait()


def parse_depths(fields: list[str], n_alt: int) -> list[tuple[int, ...] | None]:
    formats = fields[8].split(":")
    if "AD" not in formats:
        return [None] * (len(fields) - 9)
    ad_index = formats.index("AD")
    result = []
    for sample in fields[9:]:
        values = sample.split(":")
        if ad_index >= len(values) or values[ad_index] in {"", "."}:
            result.append(None)
            continue
        pieces = values[ad_index].split(",")
        if len(pieces) != n_alt + 1 or any(not value.isdigit() for value in pieces):
            result.append(None)
            continue
        counts = tuple(int(value) for value in pieces)
        result.append(counts if sum(counts) > 0 else None)
    return result


def summarize(depths: list, alt_index: int, indices: list[int]) -> dict:
    """Summarize target-ALT-positive carriers within the eligible sample cohort."""
    values, coverages = [], []
    for index in indices:
        counts = depths[index]
        if counts is None or counts[alt_index] == 0:
            continue
        depth = sum(counts)
        if depth:
            values.append(counts[alt_index] / depth)
            coverages.append(depth)
    if not values:
        return {"n_carriers": 0, "median_vaf": math.nan, "median_depth": math.nan}
    return {
        "n_carriers": len(values), "median_vaf": statistics.median(values),
        "median_depth": statistics.median(coverages),
    }


def site_failures(cohort: str, filter_value: str, stats: dict,
                  population_size: int) -> list[str]:
    """Return variant-level exclusion reasons under the cohort's fixed criteria."""
    failures = []
    if not stats["n_carriers"] or not stats["median_vaf"] < 0.4:
        failures.append("median_vaf_not_below_0.4")
    if cohort == "WGS":
        if filter_value != "PASS":
            failures.append("filter_not_PASS")
        if not stats["median_depth"] * stats["median_vaf"] >= 3:
            failures.append("median_depth_times_median_vaf_below_3")
        if stats["n_carriers"] > 36:
            failures.append("more_than_36_carriers")
    else:
        # Integer comparison preserves the strict <0.1% boundary exactly.
        if stats["n_carriers"] * 1000 >= population_size:
            failures.append("carrier_prevalence_not_below_0.001")
    return failures


def table_writer(path: Path, columns: list[str]):
    handle = path.open("w", newline="")
    writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", extrasaction="ignore")
    writer.writeheader()
    return handle, writer


def format_number(value):
    return "NA" if isinstance(value, float) and not math.isfinite(value) else value


def run(args) -> None:
    cohort = args.cohort.lower()
    output_dir = args.output_dir or Path(f"/main/results/adsp_{cohort}")
    candidates = load_candidates(args.candidates or output_dir / "candidate_variants.tsv")
    participants = load_samples(args.samples or Path(f"/main/data/adsp_{cohort}/samples.tsv"))
    analysis_samples = list(participants)
    population_size = len(analysis_samples)
    minimum_depth = {"WGS": 15, "WES": 30}[args.cohort]
    vcfs = args.vcfs
    if vcfs is None:
        directory = Path(f"/main/data/adsp_{cohort}/vcf")
        vcfs = sorted(set(directory.glob("*.vcf")) | set(directory.glob("*.vcf.gz")) |
                      set(directory.glob("*.vcf.bgz")) | set(directory.glob("*.bcf")))
    if not vcfs or len(set(path.resolve() for path in vcfs)) != len(vcfs):
        raise ValueError("Supply a nonempty, nonduplicated list of cohort VCF shards.")
    for path in vcfs:
        if not path.is_file():
            raise FileNotFoundError(f"Required VCF input is missing: {path}")
    output_dir.mkdir(parents=True, exist_ok=True)
    summary_columns = ["SNP", "SYMBOL", "Consequence", "FILTER", "n_carriers", "median_vaf", "median_depth",
                       "qc_population_size", "carrier_prevalence", "passes_qc", "exclusion_reasons"]
    observation_columns = ["sample_id", "SNP", "SYMBOL", "Consequence", "HGVSp", "ref_depth", "alt_depth",
                           "total_depth", "VAF", "CHIP_call"]
    counts = {sample: {"positive": 0, "missing": 0} for sample in analysis_samples}
    passing_variants = 0
    seen = set()
    source_samples = None
    with tempfile.TemporaryDirectory(prefix="chip_", dir=output_dir) as temporary:
        scratch = Path(temporary)
        q_handle, q_writer = table_writer(scratch / "variant_qc.tsv", summary_columns)
        o_handle, o_writer = table_writer(scratch / "allele_observations.tsv", observation_columns)
        try:
            for vcf in vcfs:
                with vcf_stream(vcf, candidates, args.reader, scratch) as (samples, lines):
                    sample_set = set(samples)
                    if source_samples is None:
                        source_samples = sample_set
                    elif source_samples != sample_set:
                        raise ValueError(f"VCF shards have different sample sets: {vcf.name}. Use a consistent cohort.")
                    missing_samples = set(analysis_samples) - sample_set
                    if missing_samples:
                        raise ValueError(f"{len(missing_samples)} analysis sample(s) are absent from {vcf.name}.")
                    index_by_sample = {sample: index for index, sample in enumerate(samples)}
                    selected_indices = [index_by_sample[sample] for sample in analysis_samples]
                    for line in lines:
                        if line.startswith("#"):
                            continue
                        fields = line.rstrip("\r\n").split("\t")
                        if len(fields) != len(samples) + 9:
                            raise ValueError(f"Malformed VCF row in {vcf.name}: sample count mismatch.")
                        chrom, pos, _, ref, alts = fields[:5]
                        alternates = alts.split(",")
                        targets = []
                        for alt_index, alt in enumerate(alternates, start=1):
                            if alt in {"*", "."} or alt.startswith("<") or "[" in alt or "]" in alt:
                                continue
                            key = canonical_variant(f"{chrom}:{pos}:{ref}:{alt}")
                            if key in candidates:
                                targets.append((key, alt_index))
                        if not targets:
                            continue
                        depths = parse_depths(fields, len(alternates))
                        if not any(depths[index] is not None for index in selected_indices):
                            raise ValueError(f"No usable FORMAT/AD values at {chrom}:{pos}; cannot evaluate CHIP candidates.")
                        for key, alt_index in targets:
                            if key in seen:
                                raise ValueError(f"Overlapping VCF records for {candidates[key]['SNP']}; supply disjoint shards.")
                            seen.add(key)
                            candidate = candidates[key]
                            stats = summarize(depths, alt_index, selected_indices)
                            failures = site_failures(args.cohort, fields[6], stats, population_size)
                            passes = not failures
                            passing_variants += passes
                            q_writer.writerow({name: format_number(value) for name, value in {
                                **candidate, "FILTER": fields[6], **stats,
                                "qc_population_size": population_size,
                                "carrier_prevalence": stats["n_carriers"] / population_size,
                                "passes_qc": int(passes), "exclusion_reasons": ";".join(failures),
                            }.items()})
                            for sample, index in zip(analysis_samples, selected_indices):
                                ad = depths[index]
                                if ad is None:
                                    counts[sample]["missing"] += passes
                                    continue
                                depth, alt_depth = sum(ad), ad[alt_index]
                                vaf = alt_depth / depth
                                call = passes and vaf > 0.02 and depth > minimum_depth
                                counts[sample]["positive"] += call
                                if alt_depth > 0:
                                    o_writer.writerow({**candidate, "sample_id": sample, "ref_depth": ad[0], "alt_depth": alt_depth,
                                                       "total_depth": depth, "VAF": vaf, "CHIP_call": int(call)})
        finally:
            q_handle.close()
            o_handle.close()
        absent = set(candidates) - seen
        if absent:
            examples = ", ".join(candidates[key]["SNP"] for key in sorted(absent)[:3])
            raise ValueError(f"{len(absent)} candidate allele(s) were absent from the supplied VCFs (e.g. {examples}). Check assembly, normalization, and shard completeness.")
        status_columns = ["sample_id", "participant_id", "CHIP", "n_qualifying_variants", "n_missing_variants"]
        s_handle, s_writer = table_writer(scratch / "status.tsv", status_columns)
        try:
            for sample in analysis_samples:
                positive = counts[sample]["positive"]
                missing = counts[sample]["missing"]
                # A qualifying positive call determines carrier status even if
                # other sites are missing. Otherwise incomplete AD gives NA.
                s_writer.writerow({"sample_id": sample, "participant_id": participants[sample],
                                   "CHIP": 1 if positive else ("NA" if missing else 0),
                                   "n_qualifying_variants": positive, "n_missing_variants": missing})
        finally:
            s_handle.close()
        settings = {
            "cohort": args.cohort, "qc_population": "unique_blood_participants",
            "qc_population_size": population_size, "candidate_count": len(candidates),
            "passing_variants": passing_variants,
            "vaf_definition": "target_ALT / sum(all_AD)",
            "variant_median_vaf_max_exclusive": 0.4,
            "variant_required_filter": "PASS" if args.cohort == "WGS" else None,
            "variant_median_depth_times_median_vaf_min_inclusive": 3 if args.cohort == "WGS" else None,
            "variant_carrier_count_max_inclusive": 36 if args.cohort == "WGS" else None,
            "variant_carrier_prevalence_max_exclusive": 0.001 if args.cohort == "WES" else None,
            "individual_vaf_min_exclusive": 0.02,
            "individual_total_allele_depth_min_exclusive": minimum_depth,
            "status_missing_rule": "NA when no positive call and any passing variant has unusable AD",
        }
        (scratch / "run_settings.json").write_text(json.dumps(settings, indent=2) + "\n")
        for name in ("variant_qc.tsv", "allele_observations.tsv", "status.tsv", "run_settings.json"):
            (scratch / name).replace(output_dir / name)
    print(f"Processed {len(candidates)} candidate alleles and {len(analysis_samples)} analysis samples.")
    print(f"Passing sites: {passing_variants}.")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cohort", required=True, choices=("WGS", "WES"))
    parser.add_argument("--candidates", type=Path)
    parser.add_argument("--vcfs", type=Path, nargs="+")
    parser.add_argument("--samples", type=Path, help="TSV of all eligible unique blood participants; default /main/data/adsp_{cohort}/samples.tsv")
    parser.add_argument("--reader", choices=("auto", "bcftools", "scan"), default="auto")
    parser.add_argument("--output-dir", type=Path)
    args = parser.parse_args()
    try:
        run(args)
    except (ValueError, OSError, ImportError, subprocess.CalledProcessError) as error:
        parser.exit(2, f"ERROR: {error}\n")


if __name__ == "__main__":
    main()
