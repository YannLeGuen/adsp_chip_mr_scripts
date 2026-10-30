#!/usr/bin/env python3
"""Synthetic checks of scientific thresholds, identity matching, and data flow.

Run: python3 /main/tests/smoke_test.py
Requires Python, numpy, and pandas. All fixtures are invented and temporary.
Exercises command-line variant selection, CHIP status, and GWAS preparation,
plus exact threshold and phenotype-input checks.
"""

from __future__ import annotations

import os
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

sys.dont_write_bytecode = True
ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT / "scripts"), str(ROOT / "scripts/chip"), str(ROOT / "scripts/mr")]

import pandas as pd
from associate_chip_ad import prepare_cohort
from call_chip_status import load_samples, site_failures
from prepare_mr_inputs import normalize, scan_gwas
from select_chip_variants import load_reference, select_variants


class ScientificSmokeTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.temporary = tempfile.TemporaryDirectory(prefix="chip_ad_smoke_")
        cls.work = Path(cls.temporary.name)
        cls.reference = pd.DataFrame([
            {"Gene name": "TET2", "Variants_category": "Frameshift, nonsense, splice", "Exact_Missense": "", "Limit": ""},
            {"Gene name": "DNMT3A", "Variants_category": "", "Exact_Missense": "R882H", "Limit": ""},
            {"Gene name": "MLL", "Variants_category": "Frameshift", "Exact_Missense": "", "Limit": "p.10-20"},
        ])
        cls.reference.to_csv(cls.work / "reference.tsv", sep="\t", index=False)
        pd.DataFrame([
            ["chr1:100:A:C", "TET2", "stop_gained", "p.Gln10Ter", "Q/*"],
            ["1:200:C:T", "DNMT3A", "missense_variant", "ENSP_TEST:p.Arg882His", "R/H"],
            ["1:300:A:G", "TET2", "stop_gained", "p.Glu100Ter", "E/*"],
        ], columns=["SNP", "SYMBOL", "Consequence", "HGVSp", "Amino_acids"]).to_csv(
            cls.work / "annotations.tsv", sep="\t", index=False)
        cls.samples = pd.DataFrame({
            "sample_id": ["s4", "s2", "s1", "s3"], "participant_id": ["p4", "p2", "p1", "p3"],
            "body_site": ["Blood"] * 4, "technical_replicate": ["0"] * 4,
        })
        cls.samples.to_csv(cls.work / "samples.tsv", sep="\t", index=False)
        header = (
            '##fileformat=VCFv4.2\n##contig=<ID=chr1,length=10000>\n'
            '##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths">\n'
            '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1\ts2\ts3\ts4\texcluded\n'
        )
        (cls.work / "cohort.vcf").write_text(header +
            'chr1\t100\t.\tA\tC\t50\tPASS\t.\tAD\t27,3\t60,0\t.\t60,0\t10,90\n'
            'chr1\t200\t.\tC\tT\t50\tPASS\t.\tAD\t30,0\t90,10\t30,0\t30,0\t100,0\n'
            'chr1\t300\t.\tA\tC,G\t50\tPASS\t.\tAD\t80,10,10\t100,0,0\t100,0,0\t100,0,0\t100,0,0\n')
        cls.cli("chip/select_chip_variants.py", "--cohort", "WGS", "--reference", cls.work / "reference.tsv",
                "--annotations", cls.work / "annotations.tsv", "--output-dir", cls.work / "chip")
        cls.cli("chip/call_chip_status.py", *cls.call_arguments())
        cls.status = pd.read_csv(cls.work / "chip/status.tsv", sep="\t", dtype=str, keep_default_na=False)
        cls.observations = pd.read_csv(cls.work / "chip/allele_observations.tsv", sep="\t")
        cls.chip_gwas = pd.DataFrame([
            ["chr1:100:A:C", .2, .02, 1e-10, .2],
            ["1:200:C:T", .01, .05, .8, .3],  # Null CHIP association, strong AD instrument.
            ["1:300:A:G", .3, .03, 1e-12, .1],
        ], columns=["SNP", "Beta_unadj", "SE_unadj", "P_unadj", "AAF"])
        cls.ad_gwas = pd.DataFrame([
            ["1", 100, "C", "A", .02, .05, .4, .2],
            ["1", 200, "T", "C", .3, .03, 1e-12, .3],
            ["1", 300, "T", "C", .02, .05, .4, .1],  # Same position, different alleles.
        ], columns=["chromosome", "base_pair_location", "effect_allele", "other_allele", "beta",
                    "standard_error", "p_value", "effect_allele_frequency"])
        cls.chip_gwas.to_csv(cls.work / "chip_gwas.csv", index=False)
        cls.ad_gwas.to_csv(cls.work / "ad_gwas.tsv", sep="\t", index=False)
        cls.cli("mr/prepare_mr_inputs.py", "--chip-gwas", cls.work / "chip_gwas.csv", "--ad-gwas",
                cls.work / "ad_gwas.tsv", "--chip-genome-wide", "--direction", "both", "--chunksize", "1",
                "--output-dir", cls.work / "mr")

    @classmethod
    def tearDownClass(cls):
        cls.temporary.cleanup()

    @classmethod
    def cli(cls, script, *args, success=True):
        result = subprocess.run([sys.executable, str(ROOT / "scripts" / script), *map(str, args)],
                                capture_output=True, text=True,
                                env={**os.environ, "PYTHONDONTWRITEBYTECODE": "1"})
        if success and result.returncode:
            raise AssertionError(result.stderr or result.stdout)
        if not success and not result.returncode:
            raise AssertionError("Expected invalid input to fail.")
        return result

    @classmethod
    def call_arguments(cls):
        return ["--cohort", "WGS", "--candidates", cls.work / "chip/candidate_variants.tsv", "--vcfs",
                cls.work / "cohort.vcf", "--samples", cls.work / "samples.tsv", "--output-dir", cls.work / "chip",
                "--reader", "scan"]

    def phenotypes(self):
        return pd.DataFrame({
            "sample_id": ["s1", "s2", "s3", "s4"], "participant_id": ["p1", "p2", "p3", "p4"],
            "AD": ["1", "0", "1", "0"], "sex": ["1", "0", "0", "1"],
            "age": ["90+", "75", "80", "65"], "APOE": ["23", "34", "33", "44"],
            "PC1": ["0.1"] * 4, "PC2": ["0.2"] * 4, "PC3": ["0.3"] * 4, "PC4": ["0.4"] * 4,
        })

    def test_chip_pipeline_matches_sample_and_participant_ids(self):
        result = self.status.set_index("sample_id")
        self.assertEqual(self.status.sample_id.tolist(), ["s4", "s2", "s1", "s3"])
        self.assertEqual(result.loc["s1", "CHIP"], "1")
        self.assertEqual(result.loc["s2", "CHIP"], "1")
        self.assertEqual(result.loc["s4", "CHIP"], "0")
        self.assertEqual(result.loc["s2", "participant_id"], "p2")

    def test_multiallelic_reads_use_target_alt_and_total_depth(self):
        row = self.observations.loc[self.observations.SNP.eq("1:300:A:G")].iloc[0]
        self.assertEqual((row.alt_depth, row.total_depth), (10, 100))
        self.assertAlmostEqual(row.VAF, .1)

    def test_missing_genotype_does_not_become_negative(self):
        row = self.status.set_index("sample_id").loc["s3"]
        self.assertEqual(row.CHIP, "NA")
        self.assertEqual(row.n_missing_variants, "1")
        self.assertEqual(row.n_qualifying_variants, "0")

    def test_site_qc_uses_selected_blood_samples_only(self):
        qc = pd.read_csv(self.work / "chip/variant_qc.tsv", sep="\t").set_index("SNP")
        row = qc.loc["1:100:A:C"]
        # The excluded fifth VCF sample has 90% VAF and must not change this site.
        self.assertEqual(row.n_carriers, 1)
        self.assertAlmostEqual(row.median_vaf, .1)

    def test_catalog_protein_range_alias_and_exact_alleles(self):
        annotations = pd.DataFrame([
            ["1:400:A:AT", "KMT2A", "frameshift_variant", "p.Arg15fs", "R/X"],
            ["1:401:A:AT", "KMT2A", "frameshift_variant", "p.Arg25fs", "R/X"],
            ["1:200:C:T", "DNMT3A", "missense_variant", "p.Arg882His", "R/H"],
            ["1:200:C:G", "DNMT3A", "missense_variant", "p.Arg882Gly", "R/G"],
        ], columns=["SNP", "SYMBOL", "Consequence", "HGVSp", "Amino_acids"])
        selected = select_variants(annotations, load_reference(self.work / "reference.tsv"))
        self.assertEqual(set(selected.SNP), {"1:400:A:AT", "1:200:C:T"})

    def test_ambiguous_gene_alias_is_rejected(self):
        reference = self.reference.copy()
        reference.loc[2, "Gene name"] = "MLL2"
        path = self.work / "ambiguous_reference.tsv"
        reference.to_csv(path, sep="\t", index=False)
        with self.assertRaisesRegex(ValueError, "MLL2.*KMT2B.*KMT2D"):
            load_reference(path)
        annotations = pd.DataFrame([
            ["1:400:A:AT", "MLL2", "frameshift_variant", "p.Arg15fs", "R/X"],
        ], columns=["SNP", "SYMBOL", "Consequence", "HGVSp", "Amino_acids"])
        with self.assertRaisesRegex(ValueError, "MLL2.*KMT2B.*KMT2D"):
            select_variants(annotations, load_reference(self.work / "reference.tsv"))
        for approved_symbol in ["KMT2B", "KMT2D"]:
            with self.subTest(approved_symbol=approved_symbol):
                reference.loc[2, "Gene name"] = approved_symbol
                reference.to_csv(path, sep="\t", index=False)
                annotations.loc[0, "SYMBOL"] = approved_symbol
                selected = select_variants(annotations, load_reference(path))
                self.assertEqual(selected.SYMBOL.tolist(), [approved_symbol])

    def test_wgs_threshold_boundaries(self):
        baseline = {"n_carriers": 36, "median_vaf": .1, "median_depth": 30}
        self.assertEqual(site_failures("WGS", "PASS", baseline, 1000), [])
        for change, reason in [({"n_carriers": 37}, "more_than_36_carriers"),
                               ({"median_depth": 29}, "median_depth_times_median_vaf_below_3"),
                               ({"median_vaf": .4}, "median_vaf_not_below_0.4")]:
            with self.subTest(change=change):
                changed = {**baseline, **change}
                self.assertIn(reason, site_failures("WGS", "PASS", changed, 1000))

    def test_wes_prevalence_is_strictly_below_one_per_thousand(self):
        stats = {"n_carriers": 2, "median_vaf": .1, "median_depth": 60}
        self.assertIn("carrier_prevalence_not_below_0.001", site_failures("WES", "PASS", stats, 2000))
        self.assertEqual(site_failures("WES", "PASS", stats, 2001), [])

    def threshold_fixture(self, cohort, depth_pairs, n_samples):
        directory = self.work / f"{cohort.lower()}_thresholds"
        directory.mkdir(exist_ok=True)
        samples = [f"threshold_{i}" for i in range(n_samples)]
        pd.DataFrame({"sample_id": samples, "participant_id": [f"participant_{i}" for i in range(n_samples)],
                      "body_site": ["Blood"] * n_samples, "technical_replicate": ["0"] * n_samples}).to_csv(
            directory / "samples.tsv", sep="\t", index=False)
        positions = list(range(1000, 1000 + len(depth_pairs)))
        pd.DataFrame({"SNP": [f"1:{position}:A:C" for position in positions],
                      "CHROM": ["1"] * len(positions), "POS": positions, "REF": ["A"] * len(positions),
                      "ALT": ["C"] * len(positions), "SYMBOL": ["TET2"] * len(positions),
                      "Consequence": ["stop_gained"] * len(positions), "HGVSp": ["p.Gln10Ter"] * len(positions),
                      "driver_rule": ["nonsense"] * len(positions)}).to_csv(
            directory / "candidates.tsv", sep="\t", index=False)
        lines = ['##fileformat=VCFv4.2', '##contig=<ID=1,length=10000>',
                 '##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths">',
                 '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t' + '\t'.join(samples)]
        for index, (position, pair) in enumerate(zip(positions, depth_pairs)):
            depths = ["100,0"] * n_samples
            depths[index] = f"{pair[0]},{pair[1]}"
            lines.append(f"1\t{position}\t.\tA\tC\t50\tPASS\t.\tAD\t" + '\t'.join(depths))
        (directory / "cohort.vcf").write_text('\n'.join(lines) + '\n')
        self.cli("chip/call_chip_status.py", "--cohort", cohort, "--candidates", directory / "candidates.tsv",
                 "--vcfs", directory / "cohort.vcf", "--samples", directory / "samples.tsv",
                 "--reader", "scan", "--output-dir", directory)
        return pd.read_csv(directory / "status.tsv", sep="\t").set_index("sample_id")

    def test_wgs_individual_depth_and_vaf_boundaries(self):
        result = self.threshold_fixture("WGS", [(12, 3), (13, 3), (147, 3), (97, 3)], 4)
        # Depth must exceed 15 and VAF must exceed 2%, both strictly.
        self.assertEqual(result.CHIP.tolist(), [0, 1, 0, 1])

    def test_wes_individual_depth_and_vaf_boundaries(self):
        result = self.threshold_fixture("WES", [(29, 1), (30, 1), (98, 2), (97, 3)], 1001)
        # Each site has one carrier in 1,001 samples, passing the prevalence rule.
        # Depth 30 and VAF exactly 2% fail; depth 31 and VAF above 2% pass.
        self.assertEqual(result.CHIP.iloc[:4].tolist(), [0, 1, 0, 1])
        self.assertEqual(int(result.CHIP.sum()), 2)

    def test_sample_manifest_requires_blood_and_no_replicates(self):
        for field, value in [("body_site", "Brain"), ("technical_replicate", "1")]:
            with self.subTest(field=field):
                samples = self.samples.copy()
                samples.loc[0, field] = value
                samples.to_csv(self.work / "invalid_samples.tsv", sep="\t", index=False)
                with self.assertRaises(ValueError):
                    load_samples(self.work / "invalid_samples.tsv")
        self.samples[["sample_id"]].to_csv(self.work / "invalid_samples.tsv", sep="\t", index=False)
        with self.assertRaises(ValueError):
            load_samples(self.work / "invalid_samples.tsv")

    def test_duplicate_participant_is_rejected_before_calling(self):
        samples = self.samples.copy()
        samples.loc[1, "participant_id"] = samples.loc[0, "participant_id"]
        samples.to_csv(self.work / "duplicates.tsv", sep="\t", index=False)
        with self.assertRaisesRegex(ValueError, "participant"):
            load_samples(self.work / "duplicates.tsv")

    def test_absent_sample_is_a_cli_error(self):
        pd.DataFrame({"sample_id": ["absent"], "participant_id": ["p_absent"],
                      "body_site": ["Blood"], "technical_replicate": ["0"]}).to_csv(
            self.work / "absent.tsv", sep="\t", index=False)
        arguments = self.call_arguments()
        arguments[arguments.index("--samples") + 1] = self.work / "absent.tsv"
        result = self.cli("chip/call_chip_status.py", *arguments, success=False)
        self.assertIn("absent from", result.stderr)

    def test_mr_normalization_preserves_effect_orientation(self):
        standardized, counts = normalize(self.ad_gwas, "AD")
        row = standardized.loc[standardized.pos.eq(200)].iloc[0]
        self.assertEqual((row.SNP, row.effect_allele, row.other_allele), ("1:200:c:t", "T", "C"))
        self.assertEqual(row.beta, .3)
        self.assertEqual(counts["invalid_rows_removed"], 0)

    def test_mr_duplicate_conflicts_are_not_silently_selected(self):
        duplicate = pd.concat([self.chip_gwas.iloc[[0]], self.chip_gwas.iloc[[0]]], ignore_index=True)
        duplicate.to_csv(self.work / "duplicate_gwas.csv", index=False)
        rows, counts = scan_gwas(self.work / "duplicate_gwas.csv", "CHIP", cutoff=5e-8, chunksize=1)
        self.assertEqual((len(rows), counts["identical_selected_duplicates_removed"]), (1, 1))
        duplicate.loc[1, "Beta_unadj"] = -.2
        duplicate.to_csv(self.work / "duplicate_gwas.csv", index=False)
        with self.assertRaisesRegex(ValueError, "conflicting summary rows"):
            scan_gwas(self.work / "duplicate_gwas.csv", "CHIP", cutoff=5e-8, chunksize=1)

    def test_mr_all_invalid_chunks_have_valid_empty_schema(self):
        for snp in ["", None, "malformed", "1:100:A:C:T"]:
            with self.subTest(snp=snp):
                invalid = self.chip_gwas.iloc[[0]].copy()
                invalid.loc[0, "SNP"] = snp
                standardized, counts = normalize(invalid, "CHIP")
                self.assertTrue(standardized.empty)
                self.assertEqual(counts["invalid_rows_removed"], 1)
                self.assertIn("effect_allele", standardized.columns)

    def test_bidirectional_mr_uses_independent_exposure_selection(self):
        forward = pd.read_csv(self.work / "mr/chip_to_ad_outcome.tsv", sep="\t")
        reverse_exp = pd.read_csv(self.work / "mr/ad_to_chip_exposure.tsv", sep="\t")
        reverse_out = pd.read_csv(self.work / "mr/ad_to_chip_outcome.tsv", sep="\t")
        self.assertEqual(forward.SNP.tolist(), ["1:100:a:c"])
        self.assertEqual(reverse_exp.SNP.tolist(), ["1:200:c:t"])
        self.assertEqual(reverse_out.SNP.tolist(), ["1:200:c:t"])
        self.assertEqual(reverse_out["pval.outcome"].item(), .8)

    def test_reverse_mr_rejects_significance_filtered_chip_input(self):
        self.chip_gwas.loc[self.chip_gwas.P_unadj.lt(5e-8)].to_csv(self.work / "filtered_gwas.csv", index=False)
        result = self.cli("mr/prepare_mr_inputs.py", "--chip-gwas", self.work / "filtered_gwas.csv", "--ad-gwas",
                          self.work / "ad_gwas.tsv", "--chip-genome-wide", "--direction", "ad_to_chip",
                          "--output-dir", self.work / "invalid_mr", success=False)
        self.assertIn("significance-filtered", result.stderr)

    def test_association_excludes_unknown_status_and_maps_apoe_dosage(self):
        frame, report = prepare_cohort(self.phenotypes(), self.status, "WGS")
        self.assertEqual((report["complete_cases"], report["missing_CHIP_status"]), (3, 1))
        row = frame.set_index("sample_id").loc["s1"]
        self.assertEqual((row.age, row.APOE2, row.APOE4), (90, 1, 0))

    def test_association_rejects_duplicate_identifiers(self):
        phenotypes = self.phenotypes()
        phenotypes.loc[1, "participant_id"] = "p1"
        with self.assertRaisesRegex(ValueError, "participant"):
            prepare_cohort(phenotypes, self.status, "WGS")
        duplicate_status = pd.concat([self.status, self.status.iloc[[0]]], ignore_index=True)
        with self.assertRaisesRegex(ValueError, "unique sample_id"):
            prepare_cohort(self.phenotypes(), duplicate_status, "WGS")

    def test_association_preserves_the_broader_qc_population(self):
        phenotypes = self.phenotypes().iloc[:3].copy()
        frame, report = prepare_cohort(phenotypes, self.status, "WGS")
        self.assertEqual(set(frame.sample_id), {"s1", "s2"})
        self.assertEqual(report["status_rows"], 4)
        self.assertEqual(report["status_samples_without_phenotype"], 1)
        self.assertEqual(report["phenotype_samples_outside_status"], 0)

    def test_association_rejects_conflicting_participant_mapping(self):
        phenotypes = self.phenotypes()
        phenotypes.loc[0, "participant_id"] = "another_participant"
        with self.assertRaisesRegex(ValueError, "does not match"):
            prepare_cohort(phenotypes, self.status, "WGS")

    def test_wes_association_restricts_european_samples(self):
        phenotypes = self.phenotypes().drop(columns=["PC3", "PC4"])
        phenotypes["ancestry"] = ["EUR", "AFR", "EUR", "EUR"]
        frame, report = prepare_cohort(phenotypes, self.status, "WES")
        self.assertEqual(set(frame.sample_id), {"s1", "s4"})
        self.assertEqual((report["ancestry_eligible"], report["complete_cases"]), (3, 2))


if __name__ == "__main__":
    unittest.main(verbosity=2)
