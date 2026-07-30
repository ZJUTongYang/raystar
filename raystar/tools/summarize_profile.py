#!/usr/bin/env python3
"""Summarize raystar_profile CSV files with nearest-rank percentiles."""

import argparse
import csv
import math
import sys
from collections import defaultdict


SCENARIO_ORDER = {
    "open_256": 0,
    "single_obstacle_256": 1,
    "narrow_gate_256": 2,
    "dense_lattice_192": 3,
    "large_open_1024": 4,
    "bundled_testmap_50": 5,
}

REQUIRED_FIELDS = {
    "schema_version",
    "scenario",
    "k",
    "phase",
    "iteration",
    "found_paths",
    "expected_paths",
    "expanded_nodes",
    "map_time_ms",
    "plan_time_ms",
    "wall_time_ms",
    "path_validation_pass",
    "terminal_consistency_pass",
    "scenario_contract_pass",
    "determinism_pass",
    "rss_available",
    "process_hwm_kib_after_plan",
    "verdict",
    "acceptance",
}


def nearest_rank(values, percentile):
    """Return the nearest-rank percentile (rank = ceil(p * n))."""

    ordered = sorted(values)
    rank = max(1, math.ceil(percentile * len(ordered)))
    return ordered[rank - 1]


def load_rows(paths):
    rows = []
    for path in paths:
        with open(path, newline="", encoding="utf-8") as stream:
            reader = csv.DictReader(stream)
            missing = REQUIRED_FIELDS.difference(reader.fieldnames or ())
            if missing:
                names = ", ".join(sorted(missing))
                raise ValueError(f"{path}: missing required CSV fields: {names}")
            for line_number, row in enumerate(reader, start=2):
                if row["schema_version"] != "2":
                    raise ValueError(
                        f"{path}:{line_number}: unsupported schema_version "
                        f"{row['schema_version']!r}"
                    )
                row["_source"] = f"{path}:{line_number}"
                row["_path"] = path
                rows.append(row)
    if not rows:
        raise ValueError("no profiling records were found")
    return rows


def summarize(rows, include_process_hwm):
    grouped = defaultdict(list)
    cases_by_path = defaultdict(set)
    for row in rows:
        try:
            key = (row["scenario"], int(row["k"]))
        except ValueError as error:
            raise ValueError(f"{row['_source']}: invalid K") from error
        grouped[key].append(row)
        cases_by_path[row["_path"]].add(key)

    if include_process_hwm:
        mixed_paths = [path for path, keys in cases_by_path.items() if len(keys) != 1]
        if mixed_paths:
            raise ValueError(
                "--process-isolated-memory requires one scenario/K per input file; "
                f"mixed inputs: {', '.join(sorted(mixed_paths))}"
            )

    summaries = []
    for (scenario, k), records in grouped.items():
        first = [row for row in records if row["phase"] == "first"]
        measured = [row for row in records if row["phase"] == "measured"]
        unknown_phases = [
            row["phase"]
            for row in records
            if row["phase"] not in {"first", "measured"}
        ]
        if len(first) != 1:
            raise ValueError(
                f"{scenario} K={k}: expected exactly one first record, got {len(first)}"
            )
        if not measured:
            raise ValueError(f"{scenario} K={k}: no measured records")
        if unknown_phases:
            raise ValueError(
                f"{scenario} K={k}: unknown phases: {sorted(set(unknown_phases))}"
            )
        try:
            first_iteration = int(first[0]["iteration"])
            measured_iterations = sorted(int(row["iteration"]) for row in measured)
        except ValueError as error:
            raise ValueError(f"{scenario} K={k}: invalid iteration") from error
        if first_iteration != 0:
            raise ValueError(f"{scenario} K={k}: first iteration must be zero")
        if measured_iterations != list(range(1, len(measured) + 1)):
            raise ValueError(
                f"{scenario} K={k}: measured iterations are not contiguous from one"
            )

        def floats(field):
            try:
                values = [float(row[field]) for row in measured]
            except ValueError as error:
                raise ValueError(f"{scenario} K={k}: invalid {field}") from error
            if any(not math.isfinite(value) or value < 0.0 for value in values):
                raise ValueError(f"{scenario} K={k}: non-finite or negative {field}")
            return values

        found = sorted({int(row["found_paths"]) for row in records})
        expected = sorted({int(row["expected_paths"]) for row in records})
        expanded = sorted({int(row["expanded_nodes"]) for row in records})
        try:
            first_wall_ms = float(first[0]["wall_time_ms"])
            hwm_values = [int(row["process_hwm_kib_after_plan"]) for row in records]
        except ValueError as error:
            raise ValueError(f"{scenario} K={k}: invalid first timing or process HWM") from error
        if not math.isfinite(first_wall_ms) or first_wall_ms < 0.0:
            raise ValueError(f"{scenario} K={k}: invalid first wall timing")
        if any(value < 0 for value in hwm_values):
            raise ValueError(f"{scenario} K={k}: negative process HWM")
        hwm_mib = max(hwm_values) / 1024.0 if include_process_hwm else None
        verdicts = sorted({row["verdict"] for row in records})
        boolean_fields = (
            "path_validation_pass",
            "terminal_consistency_pass",
            "scenario_contract_pass",
            "determinism_pass",
            "rss_available",
        )
        accepted = (
            all(row["acceptance"] == "PASS" for row in records)
            and all(row[field] == "true" for row in records for field in boolean_fields)
            and len(found) == 1
            and len(expected) == 1
            and len(expanded) == 1
        )

        map_times = floats("map_time_ms")
        plan_times = floats("plan_time_ms")
        wall_times = floats("wall_time_ms")
        summaries.append(
            {
                "scenario": scenario,
                "k": k,
                "found_paths": str(found[0]) if len(found) == 1 else f"{found[0]}..{found[-1]}",
                "expected_paths": (
                    str(expected[0])
                    if len(expected) == 1
                    else f"{expected[0]}..{expected[-1]}"
                ),
                "first_wall_ms": first_wall_ms,
                "map_p50_ms": nearest_rank(map_times, 0.50),
                "map_p95_ms": nearest_rank(map_times, 0.95),
                "plan_p50_ms": nearest_rank(plan_times, 0.50),
                "plan_p95_ms": nearest_rank(plan_times, 0.95),
                "wall_p50_ms": nearest_rank(wall_times, 0.50),
                "wall_p95_ms": nearest_rank(wall_times, 0.95),
                "expanded_nodes": (
                    str(expanded[0])
                    if len(expanded) == 1
                    else f"{expanded[0]}..{expanded[-1]}"
                ),
                "process_hwm_mib": hwm_mib,
                "measured_samples": len(measured),
                "verdict": ";".join(verdicts),
                "acceptance": "PASS" if accepted else "FAIL",
            }
        )

    return sorted(
        summaries,
        key=lambda row: (SCENARIO_ORDER.get(row["scenario"], len(SCENARIO_ORDER)), row["k"]),
    )


def write_csv(summaries):
    fields = [
        "scenario",
        "k",
        "found_paths",
        "expected_paths",
        "first_wall_ms",
        "map_p50_ms",
        "map_p95_ms",
        "plan_p50_ms",
        "plan_p95_ms",
        "wall_p50_ms",
        "wall_p95_ms",
        "expanded_nodes",
        "process_hwm_mib",
        "measured_samples",
        "verdict",
        "acceptance",
    ]
    writer = csv.DictWriter(sys.stdout, fieldnames=fields, lineterminator="\n")
    writer.writeheader()
    for summary in summaries:
        row = dict(summary)
        for field in (
            "first_wall_ms",
            "map_p50_ms",
            "map_p95_ms",
            "plan_p50_ms",
            "plan_p95_ms",
            "wall_p50_ms",
            "wall_p95_ms",
        ):
            row[field] = f"{row[field]:.3f}"
        row["process_hwm_mib"] = (
            f"{row['process_hwm_mib']:.3f}"
            if row["process_hwm_mib"] is not None
            else ""
        )
        writer.writerow(row)


def write_markdown(summaries):
    print(
        "| Scenario | K | Found / expected | First wall ms | "
        "Wall p50 / p95 ms | Map p50 / p95 ms | Plan p50 / p95 ms | "
        "Expanded nodes | Process HWM MiB | Verdict | Acceptance |"
    )
    print("|---|---:|---:|---:|---:|---:|---:|---:|---:|---|---|")
    for row in summaries:
        process_hwm = (
            f"{row['process_hwm_mib']:.3f}"
            if row["process_hwm_mib"] is not None
            else "n/a"
        )
        print(
            f"| `{row['scenario']}` | {row['k']} | "
            f"{row['found_paths']} / {row['expected_paths']} | "
            f"{row['first_wall_ms']:.3f} | "
            f"{row['wall_p50_ms']:.3f} / {row['wall_p95_ms']:.3f} | "
            f"{row['map_p50_ms']:.3f} / {row['map_p95_ms']:.3f} | "
            f"{row['plan_p50_ms']:.3f} / {row['plan_p95_ms']:.3f} | "
            f"{row['expanded_nodes']} | {process_hwm} | {row['verdict']} | "
            f"**{row['acceptance']}** |"
        )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("csv", nargs="+", help="raystar_profile CSV file(s)")
    parser.add_argument(
        "--format",
        choices=("markdown", "csv"),
        default="markdown",
        help="summary output format (default: markdown)",
    )
    parser.add_argument(
        "--process-isolated-memory",
        action="store_true",
        help=(
            "report process HWM under the caller's declaration that every scenario/K "
            "input was captured in its own fresh process"
        ),
    )
    arguments = parser.parse_args()
    summaries = summarize(load_rows(arguments.csv), arguments.process_isolated_memory)
    if arguments.format == "csv":
        write_csv(summaries)
    else:
        write_markdown(summaries)
    return 0 if all(row["acceptance"] == "PASS" for row in summaries) else 1


if __name__ == "__main__":
    try:
        sys.exit(main())
    except (OSError, ValueError) as error:
        print(f"raystar_profile_summary: {error}", file=sys.stderr)
        sys.exit(2)
