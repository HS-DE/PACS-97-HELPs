#!/usr/bin/env python3
"""Refilter ONLY pg_matrix.tsv by current HC/S + nonblank QC metadata samples.
Preserve protein rows, annotation strings and retained raw intensity strings.
Reentrant and validates output before replacing the source.
"""
import hashlib
import json
import os
import tempfile
from pathlib import Path
from yeast_protein_qc import metadata, GROUPS

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "pg_matrix.tsv"
REPORT = ROOT / "FILTERING_REPORT.md"
ANNOTATIONS = [
    "Protein.Group", "Protein.Ids", "Protein.Names", "Genes",
    "First.Protein.Description",
]
EXPECTED = {"Internal-QC": 6, "pooled QC": 15, "HC": 27, "S": 27}


def sha256(path):
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def main():
    groups = metadata()
    if len(groups) != 75:
        raise ValueError(f"Expected exactly 75 metadata samples, got {len(groups)}")
    counts = {g: sum(x == g for x in groups.values()) for g in GROUPS}
    if counts != EXPECTED:
        raise ValueError(f"Incorrect metadata group counts: {counts}")

    original_size = SRC.stat().st_size
    original_sha = sha256(SRC)
    temporary = None

    try:
        with SRC.open("r", encoding="utf-8-sig", newline="") as incoming:
            raw_header = incoming.readline()
            if not raw_header:
                raise ValueError("Empty pg_matrix.tsv")
            newline = "\r\n" if raw_header.endswith("\r\n") else "\n"
            header = raw_header.rstrip("\r\n").split("\t")
            if header[:5] != ANNOTATIONS:
                raise ValueError(
                    f"Unexpected annotation columns: {header[:5]!r}. "
                    f"Expected {ANNOTATIONS!r}"
                )
            sample_columns = header[5:]
            named = [x for x in sample_columns if x != ""]
            blank_cols = sum(x == "" for x in sample_columns)
            indices = list(range(5)) + [
                i for i, sid in enumerate(header) if i >= 5 and sid in groups
            ]
            retained = [header[i] for i in indices[5:]]
            if len(retained) != 75 or set(retained) != set(groups):
                absent = sorted(set(groups) - set(retained))
                duplicated = sorted({x for x in retained if retained.count(x) > 1})
                raise ValueError(
                    "The uploaded matrix does not contain exactly one of each "
                    f"required sample. Missing={absent}; duplicated={duplicated}"
                )
            removed = [x for x in sample_columns if x and x not in groups]
            if len(named) != len(set(named)):
                raise ValueError("Duplicate named sample IDs detected in matrix header")
            unchanged_columns = (len(header) == 80 and indices == list(range(80)))
            max_index = max(indices)
            fd, path = tempfile.mkstemp(prefix="pg_matrix_refilter_", suffix=".tsv", dir=ROOT)
            os.close(fd)
            temporary = Path(path)
            original_proteins = 0
            with temporary.open("w", encoding="utf-8", newline="") as outgoing:
                outgoing.write("\t".join(ANNOTATIONS + retained) + newline)
                for line in incoming:
                    if not line.strip():
                        raise ValueError(f"Unexpected blank protein row at line {original_proteins+2}")
                    cells = line.rstrip("\r\n").split("\t")
                    if len(cells) <= max_index:
                        raise ValueError(
                            f"Protein row {original_proteins+2} has only {len(cells)} cells "
                            f"but needs index {max_index}"
                        )
                    outgoing.write("\t".join(cells[i] for i in indices) + newline)
                    original_proteins += 1

        if original_proteins < 1:
            raise ValueError("Matrix contains no protein records")
        with temporary.open("r", encoding="utf-8", newline="") as check:
            actual_head = check.readline().rstrip("\r\n").split("\t")
            if actual_head != ANNOTATIONS + retained or len(actual_head) != 80:
                raise AssertionError("Output header validation failed")
            output_proteins = 0
            for line in check:
                if len(line.rstrip("\r\n").split("\t")) != 80:
                    raise AssertionError(
                        f"Output row {output_proteins + 2} has incorrect column count"
                    )
                output_proteins += 1
            if output_proteins != original_proteins:
                raise AssertionError("Protein row count changed")

        before_head_cols = len(header)
        if unchanged_columns:
            temporary.unlink()
            temporary = None
            action = "already filtered (no matrix rewrite)"
            filtered_sha = original_sha
            filtered_size = original_size
        else:
            temporary.replace(SRC)
            temporary = None
            action = "filtered and rewritten"
            filtered_sha = sha256(SRC)
            filtered_size = SRC.stat().st_size

        lines = [
            "# Latest pg_matrix.tsv metadata-based refilter",
            "",
            "This run filters ONLY pg_matrix.tsv, using the existing current",
            "all_sample_metadata.xlsx and all_QC_metadata.xlsx sample IDs.",
            "No other input files or prior trial result files were modified by this script.",
            "",
            f"- Status: **{action}**.",
            f"- Protein.Group rows: **{original_proteins}** before and after.",
            f"- Original TSV columns: **{before_head_cols}** (including five protein annotations).",
            "- Retained TSV columns: **80** = five protein annotations + 75 samples.",
            f"- Original named sample columns: **{len(named)}**.",
            f"- Original unnamed sample columns: **{blank_cols}**.",
            f"- Removed metadata-unlisted named columns: **{len(removed)}**.",
            "- Retained groups: **Internal-QC (6), pooled QC (15), HC (27), S (27)**.",
            "- Blank QC and other metadata-unlisted sample columns are excluded.",
            "- Sample order follows the newly uploaded matrix, without reordering.",
            "- All protein annotations and retained intensity strings are unchanged.",
            "- No normalization, log transform, or protein-row filtering is performed.",
            f"- Uploaded matrix SHA-256: {original_sha}.",
            f"- Filtered matrix SHA-256: {filtered_sha}.",
            f"- Original size: {original_size} bytes; filtered size: {filtered_size} bytes.",
            "",
            "## Excluded named sample columns",
            "",
        ]
        lines.extend([f"- {x}" for x in removed] or ["None."])
        lines.extend([
            "",
            "Note: downstream Trial 1/2 results generated from an earlier matrix",
            "should be considered historical until explicitly recomputed.",
            "",
        ])
        REPORT.write_text("\n".join(lines), encoding="utf-8")
        print(json.dumps({
            "status": "PASS", "action": action, "protein_rows": original_proteins,
            "source_column_count": before_head_cols, "output_column_count": 80,
            "source_named_samples": len(named), "removed_named_samples": removed,
            "removed_empty_columns": blank_cols, "retained_samples": len(retained),
            "groups": counts, "source_sha256": original_sha,
            "filtered_sha256": filtered_sha,
        }, ensure_ascii=False))
    finally:
        if temporary is not None and temporary.exists():
            temporary.unlink()


if __name__ == "__main__":
    main()
