#!/usr/bin/env python3
"""Select Bouzid CHIP driver candidates from ADSP variant annotations.

Annotations and the input VCF must use the same reference assembly and allele
representation. Supply one aligned row per variant/transcript. Candidate
variants are identified by their complete CHROM:POS:REF:ALT keys.
"""

from __future__ import annotations

import argparse
import re
from pathlib import Path
from urllib.parse import unquote

import pandas as pd


GENE_ALIASES = {"MLL": "KMT2A", "SFRS2": "SRSF2"}
AMINO_ACIDS = dict(zip(
    ["Ala", "Arg", "Asn", "Asp", "Cys", "Gln", "Glu", "Gly", "His", "Ile",
     "Leu", "Lys", "Met", "Phe", "Pro", "Ser", "Thr", "Trp", "Tyr", "Val", "Ter"],
    list("ARNDCQEGHILKMFPSTWYV*")
))


def canonical_gene(symbol: str) -> str:
    """Normalize unambiguous aliases and require approved symbols otherwise."""
    symbol = str(symbol).strip()
    if symbol == "MLL2":
        raise ValueError("MLL2 is an ambiguous gene symbol; supply the approved HGNC symbol KMT2B or KMT2D.")
    return GENE_ALIASES.get(symbol, symbol)


def canonical_variant(value: str) -> tuple[str, int, str, str]:
    """Normalize chromosome spelling only; never guess liftover or indel alleles."""
    parts = str(value).split(":")
    if len(parts) != 4:
        raise ValueError(f"Expected CHROM:POS:REF:ALT, got {value!r}")
    chrom, raw_pos, ref, alt = parts
    chrom = re.sub(r"^chr", "", chrom, flags=re.I)
    chrom = "MT" if chrom.upper() in {"M", "MT"} else chrom.upper()
    pos = int(raw_pos)
    ref, alt = ref.upper(), alt.upper()
    if not chrom or pos < 1 or not re.fullmatch(r"[ACGTN]+", ref):
        raise ValueError(f"Invalid sequence variant {value!r}")
    if not re.fullmatch(r"[ACGTN]+", alt) or ref == alt:
        raise ValueError(f"Expected a single, nonsymbolic alternate allele: {value!r}")
    return chrom, pos, ref, alt


def protein_notation(value: str) -> str:
    """Convert HGVSp/three-letter amino acids to comparable protein notation."""
    value = unquote(str(value)).strip()
    if ":" in value:
        value = value.rsplit(":", 1)[1]
    value = re.sub(r"^p\.", "", value).strip("()")
    for full, short in AMINO_ACIDS.items():
        value = value.replace(full, short)
    return value.replace(" ", "").replace("X", "*")


def protein_details(row: dict[str, str]) -> tuple[int | None, set[str]]:
    """Recover the position and exact amino-acid change without brittle slicing."""
    hgvsp = protein_notation(row.get("HGVSp", ""))
    changes = {hgvsp} if hgvsp not in {"", "-", ".", "?"} else set()
    match = re.match(r"[A-Z*](\d+)", hgvsp)
    position = int(match.group(1)) if match else None
    raw_position = row.get("Protein_position", "")
    if position is None and re.fullmatch(r"\d+(?:-\d+)?", raw_position):
        position = int(raw_position.split("-")[0])
    aa = row.get("Amino_acids", "")
    if position is not None and re.fullmatch(r"[A-Z*]/[A-Z*]", aa):
        before, after = aa.split("/")
        changes.add(f"{before}{position}{after}")
    return position, changes


def read_table(path: Path) -> pd.DataFrame:
    if not path.is_file():
        raise FileNotFoundError(f"Required input is missing: {path}")
    if path.suffix.lower() in {".xlsx", ".xls"}:
        frame = pd.read_excel(path, dtype=str).fillna("")
    else:
        frame = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)
    frame.columns = [str(column).strip() for column in frame.columns]
    return frame


def load_reference(path: Path) -> dict[str, list[dict]]:
    table = read_table(path)
    required = {"Gene name", "Variants_category", "Exact_Missense", "Limit"}
    missing = sorted(required - set(table.columns))
    if missing:
        raise ValueError(f"Bouzid reference lacks columns: {', '.join(missing)}")
    if table.empty:
        raise ValueError("The Bouzid reference contains no rules.")
    rules: dict[str, list[dict]] = {}
    for index, row in enumerate(table.to_dict("records"), start=2):
        input_gene = row["Gene name"].strip()
        if not input_gene:
            if not any(str(value).strip() for value in row.values()):
                continue
            raise ValueError(f"Reference row {index} has no gene name.")
        gene = canonical_gene(input_gene)
        raw_limit = row["Limit"].strip()
        limit = None
        if raw_limit and raw_limit not in {"-", "."}:
            match = re.fullmatch(r"(?:p\.)?(\d+)\s*[-–]\s*(\d+)", raw_limit)
            if not match or int(match.group(1)) > int(match.group(2)):
                raise ValueError(f"Unrecognized protein-position Limit for {gene}: {raw_limit!r}")
            limit = (int(match.group(1)), int(match.group(2)))
        exact = {
            protein_notation(value) for value in re.split(r"[,;]", row["Exact_Missense"])
            if value.strip() and value.strip() not in {"-", "."}
        }
        for change in exact:
            if not re.match(r"[A-Z*]\d+", change) or re.search(r"\s|[<>]", change):
                raise ValueError(f"Unrecognized exact protein variant for {gene}: {change!r}")
        categories = row["Variants_category"].casefold()
        if not exact and not any(term in categories for term in ("frameshift", "nonsense", "splice")):
            raise ValueError(f"Reference rule for {gene} has no supported consequence or exact variant.")
        rules.setdefault(gene, []).append({
            "categories": categories,
            "exact": exact, "limit": limit,
        })
    return rules


def match_rule(row: dict[str, str], rule: dict) -> list[str]:
    """Match curated consequence classes, protein intervals, or exact changes."""
    position, changes = protein_details(row)
    raw_consequence = row["Consequence"]
    consequences = set(re.split(r"[,&]", raw_consequence))
    limit = rule["limit"]
    in_range = limit is None or (position is not None and limit[0] <= position <= limit[1])
    reasons = []
    if in_range:
        if "frameshift" in rule["categories"] and "frameshift_variant" in consequences:
            reasons.append("frameshift")
        if "nonsense" in rule["categories"] and consequences & {"stop_gained", "stop_lost", "start_lost"}:
            reasons.append("nonsense")
        if "splice" in rule["categories"] and any(
            term.startswith("splice_") and term != "splice_polypyrimidine_tract_variant"
            for term in consequences
        ):
            reasons.append("splice")
    # The exact list may contain named indels as well as missense substitutions.
    # Such named variants are selected only when their full protein change matches.
    if changes & rule["exact"]:
        reasons.append("exact_protein_variant")
    return reasons


def select_variants(annotations: pd.DataFrame, rules: dict) -> pd.DataFrame:
    required = {"SNP", "SYMBOL", "Consequence", "HGVSp", "Amino_acids"}
    missing = sorted(required - set(annotations.columns))
    if missing:
        raise ValueError(f"Annotation table lacks columns: {', '.join(missing)}")
    selected: dict[tuple, dict] = {}
    for row_number, row in enumerate(annotations.to_dict("records"), start=2):
        input_gene = row["SYMBOL"].strip()
        # Do not cross-match genes and consequences from unaligned transcript lists.
        if "|" in input_gene or "|" in row["HGVSp"] or "|" in row["Consequence"]:
            raise ValueError(f"Annotation row {row_number}: supply one aligned variant/transcript per row, not pipe-joined annotations.")
        gene = canonical_gene(input_gene)
        if gene not in rules:
            continue
        reasons: set[str] = set()
        for rule in rules[gene]:
            reasons.update(match_rule(row, rule))
        if not reasons:
            continue
        key = canonical_variant(row["SNP"])
        output = selected.setdefault(key, {
            "SNP": ":".join(map(str, key)), "CHROM": key[0], "POS": key[1],
            "REF": key[2], "ALT": key[3], "SYMBOL": set(), "Consequence": set(),
            "HGVSp": set(), "driver_rule": set(),
        })
        output["SYMBOL"].add(gene)
        output["Consequence"].add(row["Consequence"])
        if row["HGVSp"] not in {"", "-", "."}:
            output["HGVSp"].add(row["HGVSp"])
        output["driver_rule"].update(reasons)
    if not selected:
        raise ValueError("No CHIP driver candidates matched. Check annotation columns, gene names, and the reference.")
    records = []
    for key in sorted(selected, key=lambda k: (k[0], k[1], k[2], k[3])):
        row = selected[key]
        for field in ("SYMBOL", "Consequence", "HGVSp", "driver_rule"):
            row[field] = ";".join(sorted(row[field]))
        records.append(row)
    return pd.DataFrame(records)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cohort", required=True, choices=("WGS", "WES"))
    parser.add_argument("--annotations", type=Path, help="ADSP annotation TSV; default /main/data/adsp_{cohort}/variant_annotations.tsv")
    parser.add_argument("--reference", type=Path, default=Path("/main/reference/bouzid_chip_variants.xlsx"))
    parser.add_argument("--output-dir", type=Path, help="Default /main/results/adsp_{cohort}")
    args = parser.parse_args()
    cohort = args.cohort.lower()
    annotations = args.annotations or Path(f"/main/data/adsp_{cohort}/variant_annotations.tsv")
    output_dir = args.output_dir or Path(f"/main/results/adsp_{cohort}")
    try:
        selected = select_variants(read_table(annotations), load_reference(args.reference))
        output_dir.mkdir(parents=True, exist_ok=True)
        selected.to_csv(output_dir / "candidate_variants.tsv", sep="\t", index=False)
        # This is a 1-based point list, deliberately not a BED file.
        regions = selected[["CHROM", "POS"]].drop_duplicates()
        regions.to_csv(output_dir / "candidate_regions.tsv", sep="\t", index=False, header=False)
        print(f"Selected {len(selected)} distinct CHIP driver candidates.")
    except (ValueError, OSError, ImportError) as error:
        parser.exit(2, f"ERROR: {error}\n")


if __name__ == "__main__":
    main()
