#!/usr/bin/env python3
"""Summarize raystar_profile CSV files with nearest-rank percentiles."""

import argparse
from collections import defaultdict
import csv
import math
import re
import sys


SCENARIO_ORDER = {
    "open_256": 0,
    "single_obstacle_256": 1,
    "narrow_gate_256": 2,
    "dense_lattice_192": 3,
    "large_open_1024": 4,
    "bundled_testmap_50": 5,
}

STANDARD_K_VALUES = (1, 3, 10, 50)

SCENARIO_PATH_CAPS = {
    "open_256": 1,
    "single_obstacle_256": 2,
    "narrow_gate_256": 5,
    "dense_lattice_192": None,
    "large_open_1024": 1,
    "bundled_testmap_50": None,
}

CSV_FIELDS = [
    "schema_version",
    "scenario",
    "width",
    "height",
    "resolution",
    "input_occupied_cells",
    "start_x",
    "start_y",
    "goal_x",
    "goal_y",
    "k",
    "phase",
    "iteration",
    "allow_self_crossing",
    "max_nodes",
    "timeout_ms",
    "success",
    "outcome",
    "limit_reached",
    "request_satisfied",
    "search_complete",
    "found_paths",
    "expected_paths",
    "expanded_nodes",
    "map_time_ms",
    "plan_time_ms",
    "wall_time_ms",
    "total_path_points",
    "shortest_cost_cells",
    "longest_cost_cells",
    "path_digest",
    "path_validation_pass",
    "terminal_consistency_pass",
    "scenario_contract_pass",
    "determinism_pass",
    "rss_available",
    "process_hwm_kib_after_plan",
    "verdict",
    "acceptance",
]

REQUIRED_FIELDS = set(CSV_FIELDS)


def nearest_rank(values, percentile):
    """Return the nearest-rank percentile (rank = ceil(p * n))."""
    ordered = sorted(values)
    rank = max(1, math.ceil(percentile * len(ordered)))
    return ordered[rank - 1]


def load_rows(paths):
    rows = []
    for path in paths:
        file_records = 0
        with open(path, newline="", encoding="utf-8") as stream:
            reader = csv.DictReader(stream)
            if reader.fieldnames != CSV_FIELDS:
                raise ValueError(f"{path}: header does not match profiling CSV schema v2")
            missing = REQUIRED_FIELDS.difference(reader.fieldnames or ())
            if missing:
                names = ", ".join(sorted(missing))
                raise ValueError(f"{path}: missing required CSV fields: {names}")
            for line_number, row in enumerate(reader, start=2):
                if None in row:
                    raise ValueError(f"{path}:{line_number}: unexpected extra CSV fields")
                empty = [field for field in REQUIRED_FIELDS if row[field] in (None, "")]
                if empty:
                    names = ", ".join(sorted(empty))
                    raise ValueError(
                        f"{path}:{line_number}: empty required CSV fields: {names}"
                    )
                if row["schema_version"] != "2":
                    raise ValueError(
                        f"{path}:{line_number}: unsupported schema_version "
                        f"{row['schema_version']!r}"
                    )
                row["_source"] = f"{path}:{line_number}"
                row["_path"] = path
                rows.append(row)
                file_records += 1
        if file_records == 0:
            raise ValueError(f"{path}: no profiling records were found")
    if not rows:
        raise ValueError("no profiling records were found")
    return rows


def summarize(
    rows, include_process_hwm, expected_measured_samples, require_standard_matrix
):
    grouped = defaultdict(list)
    cases_by_path = defaultdict(set)
    paths_by_case = defaultdict(set)
    for row in rows:
        try:
            key = (row["scenario"], int(row["k"]))
        except ValueError as error:
            raise ValueError(f"{row['_source']}: invalid K") from error
        if not 1 <= key[1] <= 100:
            raise ValueError(f"{row['_source']}: K must be between 1 and 100")
        if key[0] not in SCENARIO_PATH_CAPS:
            raise ValueError(f"{row['_source']}: unknown scenario {key[0]!r}")
        grouped[key].append(row)
        cases_by_path[row["_path"]].add(key)
        paths_by_case[key].add(row["_path"])

    split_cases = [key for key, paths in paths_by_case.items() if len(paths) != 1]
    if split_cases:
        names = ", ".join(
            f"{scenario}/K={k}" for scenario, k in sorted(split_cases)
        )
        raise ValueError(f"each scenario/K case must come from one input file: {names}")

    if require_standard_matrix:
        expected_cases = {
            (scenario, k)
            for scenario in SCENARIO_ORDER
            for k in STANDARD_K_VALUES
        }
        actual_cases = set(grouped)
        missing_cases = sorted(expected_cases.difference(actual_cases))
        extra_cases = sorted(actual_cases.difference(expected_cases))
        if missing_cases or extra_cases:
            details = []
            if missing_cases:
                details.append(
                    "missing "
                    + ", ".join(f"{scenario}/K={k}" for scenario, k in missing_cases)
                )
            if extra_cases:
                details.append(
                    "unexpected "
                    + ", ".join(f"{scenario}/K={k}" for scenario, k in extra_cases)
                )
            raise ValueError("standard scenario matrix mismatch: " + "; ".join(details))

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
        if len(measured) != expected_measured_samples:
            raise ValueError(
                f"{scenario} K={k}: expected {expected_measured_samples} measured "
                f"records, got {len(measured)}"
            )
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
        expected_phase_order = ["first"] + ["measured"] * expected_measured_samples
        expected_iteration_order = list(range(expected_measured_samples + 1))
        try:
            actual_iteration_order = [int(row["iteration"]) for row in records]
        except ValueError as error:
            raise ValueError(f"{scenario} K={k}: invalid iteration") from error
        if [row["phase"] for row in records] != expected_phase_order:
            raise ValueError(f"{scenario} K={k}: records are not in producer phase order")
        if actual_iteration_order != expected_iteration_order:
            raise ValueError(f"{scenario} K={k}: records are not in producer iteration order")
        if any(row["k"] != str(k) for row in records):
            raise ValueError(f"{scenario} K={k}: K is not canonically encoded")
        if any(row["allow_self_crossing"] != "false" for row in records):
            raise ValueError(f"{scenario} K={k}: unexpected self-crossing policy")
        if any(
            row["max_nodes"] != "10000" or row["timeout_ms"] != "5000"
            for row in records
        ):
            raise ValueError(f"{scenario} K={k}: unexpected profiling limits")

        def floats(field, selected_records):
            try:
                values = [float(row[field]) for row in selected_records]
            except ValueError as error:
                raise ValueError(f"{scenario} K={k}: invalid {field}") from error
            if any(not math.isfinite(value) or value < 0.0 for value in values):
                raise ValueError(f"{scenario} K={k}: non-finite or negative {field}")
            return values

        all_map_times = floats("map_time_ms", records)
        all_plan_times = floats("plan_time_ms", records)
        all_wall_times = floats("wall_time_ms", records)
        if any(
            wall_time + 0.005 < map_time + plan_time
            for map_time, plan_time, wall_time in zip(
                all_map_times, all_plan_times, all_wall_times
            )
        ):
            raise ValueError(f"{scenario} K={k}: wall timing is shorter than map plus plan")

        try:
            found = sorted({int(row["found_paths"]) for row in records})
            expected = sorted({int(row["expected_paths"]) for row in records})
            expanded = sorted({int(row["expanded_nodes"]) for row in records})
            total_path_points = sorted(
                {int(row["total_path_points"]) for row in records}
            )
            first_wall_ms = float(first[0]["wall_time_ms"])
            hwm_values = [int(row["process_hwm_kib_after_plan"]) for row in records]
        except ValueError as error:
            raise ValueError(
                f"{scenario} K={k}: invalid count, first timing, or process HWM"
            ) from error
        if any(value < 0 for value in found + expected + expanded + total_path_points):
            raise ValueError(f"{scenario} K={k}: negative count")
        if not math.isfinite(first_wall_ms) or first_wall_ms < 0.0:
            raise ValueError(f"{scenario} K={k}: invalid first wall timing")
        if any(value < 0 for value in hwm_values):
            raise ValueError(f"{scenario} K={k}: negative process HWM")
        if hwm_values != sorted(hwm_values):
            raise ValueError(f"{scenario} K={k}: process HWM decreases within one case")
        invalid_digests = [
            row["_source"]
            for row in records
            if re.fullmatch(r"[0-9a-f]{16}", row["path_digest"]) is None
        ]
        if invalid_digests:
            raise ValueError(
                f"{scenario} K={k}: invalid path digest at "
                + ", ".join(invalid_digests)
            )
        hwm_mib = max(hwm_values) / 1024.0 if include_process_hwm else None
        verdicts = sorted({row["verdict"] for row in records})
        validation_boolean_fields = (
            "path_validation_pass",
            "terminal_consistency_pass",
            "scenario_contract_pass",
            "determinism_pass",
        )
        boolean_fields = validation_boolean_fields + (
            "success",
            "request_satisfied",
            "search_complete",
            "rss_available",
        )
        invalid_booleans = [
            f"{row['_source']}:{field}={row[field]!r}"
            for row in records
            for field in boolean_fields
            if row[field] not in {"true", "false"}
        ]
        if invalid_booleans:
            raise ValueError(
                f"{scenario} K={k}: invalid boolean values: "
                + ", ".join(invalid_booleans)
            )
        deterministic_fields = (
            "scenario",
            "width",
            "height",
            "resolution",
            "input_occupied_cells",
            "start_x",
            "start_y",
            "goal_x",
            "goal_y",
            "k",
            "allow_self_crossing",
            "max_nodes",
            "timeout_ms",
            "success",
            "outcome",
            "limit_reached",
            "request_satisfied",
            "search_complete",
            "found_paths",
            "expected_paths",
            "expanded_nodes",
            "total_path_points",
            "shortest_cost_cells",
            "longest_cost_cells",
            "path_digest",
        )
        deterministic_values_match = all(
            len({row[field] for row in records}) == 1
            for field in deterministic_fields
        )
        request_satisfied_is_consistent = all(
            (row["request_satisfied"] == "true")
            == (row["outcome"] == "complete" and int(row["found_paths"]) == k)
            for row in records
        )
        scenario_cap = SCENARIO_PATH_CAPS[scenario]
        scenario_expected_paths = k if scenario_cap is None else min(k, scenario_cap)
        accepted = (
            all(row["acceptance"] == "PASS" for row in records)
            and verdicts == ["PASS"]
            and all(
                row[field] == "true"
                for row in records
                for field in validation_boolean_fields
            )
            and all(row["success"] == "true" for row in records)
            and all(row["outcome"] == "complete" for row in records)
            and all(row["limit_reached"] == "none" for row in records)
            and all(row["search_complete"] == "true" for row in records)
            and all(row["rss_available"] == "true" for row in records)
            and len(found) == 1
            and len(expected) == 1
            and len(expanded) == 1
            and found[0] == expected[0]
            and expected[0] == scenario_expected_paths
            and expanded[0] > 0
            and len(total_path_points) == 1
            and total_path_points[0] >= 2 * found[0]
            and total_path_points[0] <= found[0] * (expanded[0] + 1)
            and deterministic_values_match
            and request_satisfied_is_consistent
        )

        map_times = all_map_times[1:]
        plan_times = all_plan_times[1:]
        wall_times = all_wall_times[1:]
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
        "--expected-measured-samples",
        type=int,
        required=True,
        help="required measured-record count for every scenario/K case",
    )
    parser.add_argument(
        "--require-standard-matrix",
        action="store_true",
        help="require all six standard scenarios at K=1,3,10,50 exactly once",
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
    if arguments.expected_measured_samples <= 0:
        parser.error("--expected-measured-samples must be positive")
    summaries = summarize(
        load_rows(arguments.csv),
        arguments.process_isolated_memory,
        arguments.expected_measured_samples,
        arguments.require_standard_matrix,
    )
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
