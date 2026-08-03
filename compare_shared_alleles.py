#!/usr/bin/env python3
"""Compare nucleotide substitutions at exact-coordinate BED overlaps.

Expected BED columns in both files (tab-delimited, no header):
    chromosome  start  end  consequence  REF  ALT

The script streams `bedtools intersect` output and counts each shared genomic
coordinate once, even when one of the inputs contains duplicate records or
multiple alternate alleles at that coordinate.
"""

from __future__ import annotations

import argparse
import csv
import shutil
import subprocess
import sys
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Compare REF>ALT substitutions at exact-coordinate overlaps "
            "between germline and somatic six-column BED files."
        )
    )
    parser.add_argument(
        "--germline",
        required=True,
        type=Path,
        help="Sorted six-column germline BED file.",
    )
    parser.add_argument(
        "--somatic",
        required=True,
        type=Path,
        help="Sorted six-column somatic BED file.",
    )
    parser.add_argument(
        "--prefix",
        default="germline_vs_somatic_WG",
        help="Prefix for the two output TSV files.",
    )
    return parser.parse_args()


def validate_input(path: Path) -> None:
    if not path.is_file():
        raise SystemExit(f"ERROR: file not found: {path}")
    if path.stat().st_size == 0:
        raise SystemExit(f"ERROR: file is empty: {path}")


def format_substitutions(values: set[tuple[str, str]]) -> str:
    return ";".join(f"{ref}>{alt}" for ref, alt in sorted(values))


def classify_locus(
    germline: set[tuple[str, str]], somatic: set[tuple[str, str]]
) -> tuple[str, set[tuple[str, str]], bool, bool]:
    exact = germline & somatic
    shared_ref = bool({ref for ref, _ in germline} & {ref for ref, _ in somatic})
    shared_alt = bool({alt for _, alt in germline} & {alt for _, alt in somatic})

    if exact:
        if germline == somatic:
            category = "exact_match_only"
        else:
            category = "exact_match_with_additional_substitution"
    elif shared_ref:
        category = "coordinate_only_same_REF_different_ALT"
    else:
        category = "reference_mismatch"

    return category, exact, shared_ref, shared_alt


def main() -> None:
    args = parse_args()
    validate_input(args.germline)
    validate_input(args.somatic)

    if shutil.which("bedtools") is None:
        raise SystemExit(
            "ERROR: bedtools was not found. Install it first (for example: "
            "brew install bedtools)."
        )

    detail_path = Path(f"{args.prefix}.allele_concordance_by_locus.tsv")
    summary_path = Path(f"{args.prefix}.allele_concordance_summary.tsv")

    command = [
        "bedtools",
        "intersect",
        "-sorted",
        "-wa",
        "-wb",
        "-a",
        str(args.germline),
        "-b",
        str(args.somatic),
    ]

    counts = {
        "shared_loci": 0,
        "exact_REF_ALT_match": 0,
        "no_exact_REF_ALT_match": 0,
        "exact_match_only": 0,
        "exact_match_with_additional_substitution": 0,
        "coordinate_only_same_REF_different_ALT": 0,
        "reference_mismatch": 0,
        "ALT_match": 0,
        "shared_REF": 0,
    }

    current_key: tuple[str, int, int] | None = None
    germline_substitutions: set[tuple[str, str]] = set()
    somatic_substitutions: set[tuple[str, str]] = set()

    with detail_path.open("w", newline="", encoding="utf-8") as detail_handle:
        writer = csv.writer(detail_handle, delimiter="\t", lineterminator="\n")
        writer.writerow(
            [
                "chromosome",
                "start",
                "end",
                "germline_REF_ALT",
                "somatic_REF_ALT",
                "exact_REF_ALT",
                "classification",
                "exact_match_present",
                "shared_REF_present",
                "ALT_match_present",
            ]
        )

        def write_current_locus() -> None:
            if current_key is None:
                return

            category, exact, shared_ref, shared_alt = classify_locus(
                germline_substitutions, somatic_substitutions
            )
            counts["shared_loci"] += 1
            counts[category] += 1

            if exact:
                counts["exact_REF_ALT_match"] += 1
            else:
                counts["no_exact_REF_ALT_match"] += 1
            if shared_ref:
                counts["shared_REF"] += 1
            if shared_alt:
                counts["ALT_match"] += 1

            writer.writerow(
                [
                    current_key[0],
                    current_key[1],
                    current_key[2],
                    format_substitutions(germline_substitutions),
                    format_substitutions(somatic_substitutions),
                    format_substitutions(exact) if exact else "-",
                    category,
                    int(bool(exact)),
                    int(shared_ref),
                    int(shared_alt),
                ]
            )

        process = subprocess.Popen(
            command,
            stdout=subprocess.PIPE,
            text=True,
            encoding="utf-8",
        )
        assert process.stdout is not None

        for line_number, line in enumerate(process.stdout, start=1):
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 12:
                process.kill()
                raise SystemExit(
                    "ERROR: expected 12 columns in the bedtools output "
                    f"(six per input), but found {len(fields)} at output line "
                    f"{line_number}. Confirm that both BED files contain "
                    "exactly six tab-delimited columns and no header."
                )

            g_chrom, g_start, g_end = fields[0], int(fields[1]), int(fields[2])
            s_chrom, s_start, s_end = fields[6], int(fields[7]), int(fields[8])

            # The inputs contain one-base SNVs, but explicitly enforce identical
            # coordinates so interval overlap is never mistaken for locus identity.
            if (g_chrom, g_start, g_end) != (s_chrom, s_start, s_end):
                continue

            key = (g_chrom, g_start, g_end)
            if current_key is not None and key != current_key:
                write_current_locus()
                germline_substitutions.clear()
                somatic_substitutions.clear()

            current_key = key
            germline_substitutions.add((fields[4].upper(), fields[5].upper()))
            somatic_substitutions.add((fields[10].upper(), fields[11].upper()))

        write_current_locus()

        return_code = process.wait()
        if return_code != 0:
            raise SystemExit(
                f"ERROR: bedtools intersect exited with status {return_code}."
            )

    total = counts["shared_loci"]
    if total == 0:
        detail_path.unlink(missing_ok=True)
        raise SystemExit(
            "ERROR: no exact-coordinate overlaps were found. Check genome build, "
            "chromosome naming, coordinates, and sorting."
        )

    metric_order = [
        "shared_loci",
        "exact_REF_ALT_match",
        "no_exact_REF_ALT_match",
        "exact_match_only",
        "exact_match_with_additional_substitution",
        "coordinate_only_same_REF_different_ALT",
        "reference_mismatch",
        "shared_REF",
        "ALT_match",
    ]

    with summary_path.open("w", newline="", encoding="utf-8") as summary_handle:
        writer = csv.writer(summary_handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["metric", "n_loci", "percentage_of_shared_loci"])
        for metric in metric_order:
            count = counts[metric]
            percentage = 100.0 if metric == "shared_loci" else 100.0 * count / total
            writer.writerow([metric, count, f"{percentage:.4f}"])

    concordant = counts["exact_REF_ALT_match"]
    print(f"Shared exact-coordinate loci: {total:,}")
    print(
        "Loci with at least one identical REF>ALT substitution: "
        f"{concordant:,} ({100.0 * concordant / total:.4f}%)"
    )
    print(
        "Coordinate-only loci without an identical REF>ALT substitution: "
        f"{counts['no_exact_REF_ALT_match']:,} "
        f"({100.0 * counts['no_exact_REF_ALT_match'] / total:.4f}%)"
    )
    print(f"Detailed output: {detail_path}")
    print(f"Summary output:  {summary_path}")


if __name__ == "__main__":
    main()
