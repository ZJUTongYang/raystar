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

MODE_ORDER = {
    "top_k": 0,
    "all_within_length": 1,
    "multi_goal": 2,
}

SCENARIO_PATH_CAPS = {
    "open_256": 1,
    "single_obstacle_256": 2,
    "narrow_gate_256": 5,
    "dense_lattice_192": None,
    "large_open_1024": 1,
    "bundled_testmap_50": None,
}

CSV_FIELDS_V2 = [
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

CSV_FIELDS_V3 = CSV_FIELDS_V2 + [
    "mode",
    "max_path_cost_cells",
    "goal_count",
    "completion",
    "per_goal_complete",
    "per_goal_outcomes",
    "per_goal_limits",
    "per_goal_completions",
    "per_goal_found_paths",
    "max_cost_bounded_paths",
    "max_path_points",
    "max_multi_goal_count",
]

REQUIRED_FIELDS_V2 = set(CSV_FIELDS_V2)
REQUIRED_FIELDS_V3 = set(CSV_FIELDS_V3)

VALID_OUTCOMES = {
    "invalid_request",
    "complete",
    "no_path",
    "limit_reached",
    "failed",
}

VALID_LIMITS = {
    "none",
    "max_nodes",
    "timeout",
    "cancelled",
    "max_path_points",
    "max_paths",
}

VALID_SINGLE_COMPLETIONS = {
    "none",
    "requested_k_reached",
    "frontier_exhausted",
    "cost_bound_exhausted",
}

VALID_MULTI_COMPLETIONS = {"all_goals_complete", "incomplete"}


def _legacy_v2_defaults(row):
    """Return schema-v3-equivalent fields inferred from one schema-v2 row."""
    top_k_satisfied = (
        row["outcome"] == "complete"
        and row["limit_reached"] == "none"
        and row["found_paths"] == row["k"]
    )
    complete = (
        row["outcome"] in {"complete", "no_path"}
        and row["limit_reached"] == "none"
    )
    completion = "requested_k_reached" if top_k_satisfied else "frontier_exhausted"
    if not complete:
        completion = "none"
    return {
        "mode": "top_k",
        "max_path_cost_cells": "0.000",
        "goal_count": "1",
        "completion": completion,
        "per_goal_complete": "true" if complete else "false",
        "per_goal_outcomes": row["outcome"],
        "per_goal_limits": row["limit_reached"],
        "per_goal_completions": completion,
        "per_goal_found_paths": row["found_paths"],
        "max_cost_bounded_paths": "1000",
        "max_path_points": "100000",
        "max_multi_goal_count": "32",
    }


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
            if reader.fieldnames == CSV_FIELDS_V2:
                schema_version = "2"
                required_fields = REQUIRED_FIELDS_V2
            elif reader.fieldnames == CSV_FIELDS_V3:
                schema_version = "3"
                required_fields = REQUIRED_FIELDS_V3
            else:
                raise ValueError(
                    f"{path}: header does not match profiling CSV schema v2 or v3"
                )
            missing = required_fields.difference(reader.fieldnames or ())
            if missing:
                names = ", ".join(sorted(missing))
                raise ValueError(f"{path}: missing required CSV fields: {names}")
            for line_number, row in enumerate(reader, start=2):
                if None in row:
                    raise ValueError(f"{path}:{line_number}: unexpected extra CSV fields")
                empty = [
                    field for field in required_fields if row[field] in (None, "")
                ]
                if empty:
                    names = ", ".join(sorted(empty))
                    raise ValueError(
                        f"{path}:{line_number}: empty required CSV fields: {names}"
                    )
                if row["schema_version"] != schema_version:
                    raise ValueError(
                        f"{path}:{line_number}: unsupported schema_version "
                        f"{row['schema_version']!r}"
                    )
                if schema_version == "2":
                    row.update(_legacy_v2_defaults(row))
                row["_csv_schema"] = schema_version
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
        mode = row["mode"]
        if mode not in MODE_ORDER:
            raise ValueError(f"{row['_source']}: unknown profiling mode {mode!r}")
        try:
            k = int(row["k"])
            max_path_cost = float(row["max_path_cost_cells"])
            goal_count = int(row["goal_count"])
            max_cost_bounded_paths = int(row["max_cost_bounded_paths"])
            max_path_points = int(row["max_path_points"])
            max_multi_goal_count = int(row["max_multi_goal_count"])
            max_nodes = int(row["max_nodes"])
            timeout_ms = int(row["timeout_ms"])
        except ValueError as error:
            raise ValueError(
                f"{row['_source']}: invalid mode parameter or profiling limit"
            ) from error
        if row["k"] != str(k):
            raise ValueError(f"{row['_source']}: K is not canonically encoded")
        if row["goal_count"] != str(goal_count):
            raise ValueError(
                f"{row['_source']}: goal_count is not canonically encoded"
            )
        if row["max_cost_bounded_paths"] != str(max_cost_bounded_paths):
            raise ValueError(
                f"{row['_source']}: max_cost_bounded_paths is not canonically encoded"
            )
        if row["max_path_points"] != str(max_path_points):
            raise ValueError(
                f"{row['_source']}: max_path_points is not canonically encoded"
            )
        if row["max_multi_goal_count"] != str(max_multi_goal_count):
            raise ValueError(
                f"{row['_source']}: max_multi_goal_count is not canonically encoded"
            )
        if row["max_nodes"] != str(max_nodes) or row["timeout_ms"] != str(timeout_ms):
            raise ValueError(
                f"{row['_source']}: profiling limits are not canonically encoded"
            )
        if not math.isfinite(max_path_cost) or max_path_cost < 0.0:
            raise ValueError(
                f"{row['_source']}: max_path_cost_cells must be finite and nonnegative"
            )
        if (
            max_nodes <= 0
            or timeout_ms < 0
            or max_cost_bounded_paths <= 0
        ):
            raise ValueError(f"{row['_source']}: invalid nonpositive profiling limit")
        if max_path_points < 2:
            raise ValueError(
                f"{row['_source']}: max_path_points must be at least 2"
            )
        if max_multi_goal_count <= 0:
            raise ValueError(
                f"{row['_source']}: max_multi_goal_count must be positive"
            )
        if row["_csv_schema"] == "2" and (
            max_nodes != 10000 or timeout_ms != 5000
        ):
            raise ValueError(f"{row['_source']}: unexpected schema-v2 profiling limits")
        if mode == "top_k":
            if not 1 <= k <= 100:
                raise ValueError(f"{row['_source']}: K must be between 1 and 100")
            if max_path_cost != 0.0 or goal_count != 1:
                raise ValueError(
                    f"{row['_source']}: top_k requires zero bound and goal_count=1"
                )
        elif mode == "all_within_length":
            if k != 0 or max_path_cost <= 0.0 or goal_count != 1:
                raise ValueError(
                    f"{row['_source']}: all_within_length requires K=0, a positive "
                    "bound, and goal_count=1"
                )
        elif k != 0 or max_path_cost <= 0.0 or not 2 <= goal_count <= 32:
            raise ValueError(
                f"{row['_source']}: multi_goal requires K=0, a positive bound, and "
                "goal_count between 2 and 32"
            )
        if mode == "multi_goal" and max_multi_goal_count < goal_count:
            raise ValueError(
                f"{row['_source']}: max_multi_goal_count must be at least goal_count"
            )
        scenario = row["scenario"]
        if scenario not in SCENARIO_PATH_CAPS:
            raise ValueError(f"{row['_source']}: unknown scenario {scenario!r}")
        key = (scenario, mode, k, max_path_cost, goal_count)
        grouped[key].append(row)
        cases_by_path[row["_path"]].add(key)
        paths_by_case[key].add(row["_path"])

    split_cases = [key for key, paths in paths_by_case.items() if len(paths) != 1]
    if split_cases:
        names = ", ".join(
            f"{scenario}/{mode}/K={k}/bound={bound:g}/goals={goal_count}"
            for scenario, mode, k, bound, goal_count in sorted(split_cases)
        )
        raise ValueError(f"each profiling case must come from one input file: {names}")

    if require_standard_matrix:
        expected_cases = {
            (scenario, "top_k", k, 0.0, 1)
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
                    + ", ".join(
                        f"{scenario}/{mode}/K={k}/bound={bound:g}/goals={goal_count}"
                        for scenario, mode, k, bound, goal_count in missing_cases
                    )
                )
            if extra_cases:
                details.append(
                    "unexpected "
                    + ", ".join(
                        f"{scenario}/{mode}/K={k}/bound={bound:g}/goals={goal_count}"
                        for scenario, mode, k, bound, goal_count in extra_cases
                    )
                )
            raise ValueError("standard scenario matrix mismatch: " + "; ".join(details))

    if include_process_hwm:
        mixed_paths = [path for path, keys in cases_by_path.items() if len(keys) != 1]
        if mixed_paths:
            raise ValueError(
                "--process-isolated-memory requires one profiling case per input file; "
                f"mixed inputs: {', '.join(sorted(mixed_paths))}"
            )

    summaries = []
    for (scenario, mode, k, max_path_cost, goal_count), records in grouped.items():
        case = (
            f"{scenario} mode={mode} K={k} bound={max_path_cost:g} "
            f"goals={goal_count}"
        )
        first = [row for row in records if row["phase"] == "first"]
        measured = [row for row in records if row["phase"] == "measured"]
        unknown_phases = [
            row["phase"]
            for row in records
            if row["phase"] not in {"first", "measured"}
        ]
        if len(first) != 1:
            raise ValueError(
                f"{case}: expected exactly one first record, got {len(first)}"
            )
        if len(measured) != expected_measured_samples:
            raise ValueError(
                f"{case}: expected {expected_measured_samples} measured "
                f"records, got {len(measured)}"
            )
        if unknown_phases:
            raise ValueError(
                f"{case}: unknown phases: {sorted(set(unknown_phases))}"
            )
        try:
            first_iteration = int(first[0]["iteration"])
            measured_iterations = sorted(int(row["iteration"]) for row in measured)
        except ValueError as error:
            raise ValueError(f"{case}: invalid iteration") from error
        if first_iteration != 0:
            raise ValueError(f"{case}: first iteration must be zero")
        if measured_iterations != list(range(1, len(measured) + 1)):
            raise ValueError(
                f"{case}: measured iterations are not contiguous from one"
            )
        expected_phase_order = ["first"] + ["measured"] * expected_measured_samples
        expected_iteration_order = list(range(expected_measured_samples + 1))
        try:
            actual_iteration_order = [int(row["iteration"]) for row in records]
        except ValueError as error:
            raise ValueError(f"{case}: invalid iteration") from error
        if [row["phase"] for row in records] != expected_phase_order:
            raise ValueError(f"{case}: records are not in producer phase order")
        if actual_iteration_order != expected_iteration_order:
            raise ValueError(f"{case}: records are not in producer iteration order")
        if any(row["allow_self_crossing"] != "false" for row in records):
            raise ValueError(f"{case}: unexpected self-crossing policy")

        def floats(field, selected_records):
            try:
                values = [float(row[field]) for row in selected_records]
            except ValueError as error:
                raise ValueError(f"{case}: invalid {field}") from error
            if any(not math.isfinite(value) or value < 0.0 for value in values):
                raise ValueError(f"{case}: non-finite or negative {field}")
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
            raise ValueError(f"{case}: wall timing is shorter than map plus plan")

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
                f"{case}: invalid count, first timing, or process HWM"
            ) from error
        if any(value < 0 for value in found + expected + expanded + total_path_points):
            raise ValueError(f"{case}: negative count")
        if not math.isfinite(first_wall_ms) or first_wall_ms < 0.0:
            raise ValueError(f"{case}: invalid first wall timing")
        if any(value < 0 for value in hwm_values):
            raise ValueError(f"{case}: negative process HWM")
        if hwm_values != sorted(hwm_values):
            raise ValueError(f"{case}: process HWM decreases within one case")
        invalid_digests = [
            row["_source"]
            for row in records
            if re.fullmatch(r"[0-9a-f]{16}", row["path_digest"]) is None
        ]
        if invalid_digests:
            raise ValueError(
                f"{case}: invalid path digest at "
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
                f"{case}: invalid boolean values: "
                + ", ".join(invalid_booleans)
            )

        per_goal_data = []
        for row in records:
            fields = {
                field: row[field].split("|")
                for field in (
                    "per_goal_complete",
                    "per_goal_outcomes",
                    "per_goal_limits",
                    "per_goal_completions",
                    "per_goal_found_paths",
                )
            }
            wrong_lengths = [
                field for field, values in fields.items() if len(values) != goal_count
            ]
            if wrong_lengths:
                raise ValueError(
                    f"{row['_source']}: per-goal array length does not match "
                    f"goal_count for {', '.join(wrong_lengths)}"
                )
            if any(value == "" for values in fields.values() for value in values):
                raise ValueError(f"{row['_source']}: empty per-goal array item")
            if any(
                value not in {"true", "false"}
                for value in fields["per_goal_complete"]
            ):
                raise ValueError(f"{row['_source']}: invalid per_goal_complete value")
            if any(value not in VALID_OUTCOMES for value in fields["per_goal_outcomes"]):
                raise ValueError(f"{row['_source']}: invalid per_goal_outcomes value")
            if any(value not in VALID_LIMITS for value in fields["per_goal_limits"]):
                raise ValueError(f"{row['_source']}: invalid per_goal_limits value")
            if any(
                value not in VALID_SINGLE_COMPLETIONS
                for value in fields["per_goal_completions"]
            ):
                raise ValueError(f"{row['_source']}: invalid per_goal_completions value")
            try:
                per_goal_found = [
                    int(value) for value in fields["per_goal_found_paths"]
                ]
            except ValueError as error:
                raise ValueError(
                    f"{row['_source']}: invalid per_goal_found_paths value"
                ) from error
            if any(value < 0 for value in per_goal_found):
                raise ValueError(f"{row['_source']}: negative per-goal path count")
            if any(
                encoded != str(value)
                for encoded, value in zip(
                    fields["per_goal_found_paths"], per_goal_found
                )
            ):
                raise ValueError(
                    f"{row['_source']}: per_goal_found_paths is not canonically encoded"
                )
            per_goal_data.append(
                {
                    "complete": fields["per_goal_complete"],
                    "outcomes": fields["per_goal_outcomes"],
                    "limits": fields["per_goal_limits"],
                    "completions": fields["per_goal_completions"],
                    "found": per_goal_found,
                }
            )

        invalid_terminal_enums = [
            row["_source"]
            for row in records
            if row["outcome"] not in VALID_OUTCOMES
            or row["limit_reached"] not in VALID_LIMITS
            or (
                mode == "multi_goal"
                and row["completion"] not in VALID_MULTI_COMPLETIONS
            )
            or (
                mode != "multi_goal"
                and row["completion"] not in VALID_SINGLE_COMPLETIONS
            )
        ]
        if invalid_terminal_enums:
            raise ValueError(
                f"{case}: invalid terminal enum at "
                + ", ".join(invalid_terminal_enums)
            )
        invalid_acceptance = [
            row["_source"] for row in records if row["acceptance"] not in {"PASS", "FAIL"}
        ]
        if invalid_acceptance:
            raise ValueError(
                f"{case}: invalid acceptance value at "
                + ", ".join(invalid_acceptance)
            )

        shortest_costs = floats("shortest_cost_cells", records)
        longest_costs = floats("longest_cost_cells", records)
        cost_ranges_valid = all(
            shortest <= longest
            and (
                (int(row["found_paths"]) == 0 and shortest == 0.0 and longest == 0.0)
                or (int(row["found_paths"]) > 0 and shortest > 0.0)
            )
            for row, shortest, longest in zip(
                records, shortest_costs, longest_costs
            )
        )
        bounded_costs_valid = mode == "top_k" or all(
            longest <= max_path_cost + 0.0005 for longest in longest_costs
        )
        resource_bounds_valid = all(
            int(row["expanded_nodes"]) <= int(row["max_nodes"])
            and int(row["total_path_points"]) <= int(row["max_path_points"])
            for row in records
        )

        deterministic_fields = (
            "_csv_schema",
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
            "mode",
            "max_path_cost_cells",
            "goal_count",
            "completion",
            "per_goal_complete",
            "per_goal_outcomes",
            "per_goal_limits",
            "per_goal_completions",
            "per_goal_found_paths",
            "max_cost_bounded_paths",
            "max_path_points",
            "max_multi_goal_count",
        )
        deterministic_values_match = all(
            len({row[field] for row in records}) == 1
            for field in deterministic_fields
        )
        terminal_semantics_valid = True
        per_goal_path_caps_valid = True
        for row, per_goal in zip(records, per_goal_data):
            row_found = int(row["found_paths"])
            aggregate_outcome = "complete" if row_found > 0 else "no_path"
            success_matches_paths = (row["success"] == "true") == (row_found > 0)
            if mode == "top_k":
                completion_valid = (
                    row["completion"] == "requested_k_reached"
                    and row_found == k
                ) or (
                    row["completion"] == "frontier_exhausted" and row_found < k
                )
                request_satisfied = (
                    row["completion"] == "requested_k_reached" and row_found == k
                )
                per_goal_valid = (
                    per_goal["complete"] == ["true"]
                    and per_goal["outcomes"] == [row["outcome"]]
                    and per_goal["limits"] == [row["limit_reached"]]
                    and per_goal["completions"] == [row["completion"]]
                    and per_goal["found"] == [row_found]
                )
            elif mode == "all_within_length":
                completion_valid = row["completion"] in {
                    "cost_bound_exhausted",
                    "frontier_exhausted",
                }
                request_satisfied = completion_valid
                per_goal_valid = (
                    per_goal["complete"] == ["true"]
                    and per_goal["outcomes"] == [row["outcome"]]
                    and per_goal["limits"] == [row["limit_reached"]]
                    and per_goal["completions"] == [row["completion"]]
                    and per_goal["found"] == [row_found]
                )
            else:
                per_goal_certificates = [
                    complete == "true"
                    and outcome in {"complete", "no_path"}
                    and limit == "none"
                    and completion
                    in {"cost_bound_exhausted", "frontier_exhausted"}
                    and ((found_count > 0) == (outcome == "complete"))
                    for complete, outcome, limit, completion, found_count in zip(
                        per_goal["complete"],
                        per_goal["outcomes"],
                        per_goal["limits"],
                        per_goal["completions"],
                        per_goal["found"],
                    )
                ]
                completion_valid = (
                    row["completion"] == "all_goals_complete"
                    and all(per_goal_certificates)
                )
                request_satisfied = completion_valid
                per_goal_valid = (
                    sum(per_goal["found"]) == row_found
                    and all(per_goal_certificates)
                )
            terminal_semantics_valid = terminal_semantics_valid and (
                row["outcome"] == aggregate_outcome
                and row["limit_reached"] == "none"
                and row["search_complete"] == "true"
                and success_matches_paths
                and completion_valid
                and (row["request_satisfied"] == "true") == request_satisfied
                and per_goal_valid
            )
            if mode != "top_k":
                per_goal_path_caps_valid = per_goal_path_caps_valid and all(
                    found_count <= int(row["max_cost_bounded_paths"])
                    for found_count in per_goal["found"]
                )

        scenario_cap = SCENARIO_PATH_CAPS[scenario]
        scenario_expected_paths = None
        if mode == "top_k":
            scenario_expected_paths = k if scenario_cap is None else min(k, scenario_cap)
        stable_counts = len(found) == len(expected) == len(expanded) == 1
        stable_points = len(total_path_points) == 1
        path_point_contract = stable_counts and stable_points and (
            (found[0] == 0 and total_path_points[0] == 0)
            or (
                found[0] > 0
                and total_path_points[0] >= 2 * found[0]
                and total_path_points[0] <= found[0] * (expanded[0] + 1)
            )
        )
        expected_count_valid = (
            stable_counts
            and found[0] == expected[0]
            and (
                scenario_expected_paths is None
                or expected[0] == scenario_expected_paths
            )
        )
        expanded_count_valid = stable_counts and (
            expanded[0] > 0 if mode == "top_k" else expanded[0] >= 0
        )
        accepted = (
            all(row["acceptance"] == "PASS" for row in records)
            and verdicts == ["PASS"]
            and all(
                row[field] == "true"
                for row in records
                for field in validation_boolean_fields
            )
            and all(row["rss_available"] == "true" for row in records)
            and expected_count_valid
            and expanded_count_valid
            and path_point_contract
            and cost_ranges_valid
            and bounded_costs_valid
            and resource_bounds_valid
            and per_goal_path_caps_valid
            and deterministic_values_match
            and terminal_semantics_valid
        )

        map_times = all_map_times[1:]
        plan_times = all_plan_times[1:]
        wall_times = all_wall_times[1:]
        summaries.append(
            {
                "scenario": scenario,
                "mode": mode,
                "k": k,
                "max_path_cost_cells": records[0]["max_path_cost_cells"],
                "goal_count": goal_count,
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
        key=lambda row: (
            SCENARIO_ORDER.get(row["scenario"], len(SCENARIO_ORDER)),
            MODE_ORDER[row["mode"]],
            row["k"],
            float(row["max_path_cost_cells"]),
            row["goal_count"],
        ),
    )


def write_csv(summaries):
    fields = [
        "scenario",
        "mode",
        "k",
        "max_path_cost_cells",
        "goal_count",
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
        "| Scenario | Mode | K | Bound (cells) | Goals | Found / expected | First wall ms | "
        "Wall p50 / p95 ms | Map p50 / p95 ms | Plan p50 / p95 ms | "
        "Expanded nodes | Process HWM MiB | Verdict | Acceptance |"
    )
    print(
        "|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---|---|"
    )
    for row in summaries:
        process_hwm = (
            f"{row['process_hwm_mib']:.3f}"
            if row["process_hwm_mib"] is not None
            else "n/a"
        )
        print(
            f"| `{row['scenario']}` | `{row['mode']}` | {row['k']} | "
            f"{row['max_path_cost_cells']} | {row['goal_count']} | "
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
        help="required measured-record count for every profiling case",
    )
    parser.add_argument(
        "--require-standard-matrix",
        action="store_true",
        help=(
            "require the top-k matrix of all six standard scenarios at "
            "K=1,3,10,50 exactly once"
        ),
    )
    parser.add_argument(
        "--process-isolated-memory",
        action="store_true",
        help=(
            "report process HWM under the caller's declaration that every profiling "
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
