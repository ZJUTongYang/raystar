#!/usr/bin/env python3
"""Contract tests for the installed profiling CSV summarizer."""

import csv
import pathlib
import subprocess
import sys
import tempfile
import unittest


SUMMARY = pathlib.Path(sys.argv.pop(1)).resolve()

FIELDS_V2 = [
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

FIELDS_V3 = FIELDS_V2 + [
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


def make_rows(measured_samples=2):
    base = {
        "schema_version": "2",
        "scenario": "open_256",
        "width": "256",
        "height": "256",
        "resolution": "1.000",
        "input_occupied_cells": "0",
        "start_x": "4.500",
        "start_y": "4.500",
        "goal_x": "250.500",
        "goal_y": "250.500",
        "k": "1",
        "phase": "first",
        "iteration": "0",
        "allow_self_crossing": "false",
        "max_nodes": "10000",
        "timeout_ms": "5000",
        "success": "true",
        "outcome": "complete",
        "limit_reached": "none",
        "request_satisfied": "true",
        "search_complete": "true",
        "found_paths": "1",
        "expected_paths": "1",
        "expanded_nodes": "1",
        "map_time_ms": "10.0",
        "plan_time_ms": "1.0",
        "wall_time_ms": "11.0",
        "total_path_points": "2",
        "shortest_cost_cells": "347.896",
        "longest_cost_cells": "347.896",
        "path_digest": "0123456789abcdef",
        "path_validation_pass": "true",
        "terminal_consistency_pass": "true",
        "scenario_contract_pass": "true",
        "determinism_pass": "true",
        "rss_available": "true",
        "process_hwm_kib_after_plan": "4096",
        "verdict": "PASS",
        "acceptance": "PASS",
    }
    rows = [base]
    for iteration in range(1, measured_samples + 1):
        row = dict(base)
        row["phase"] = "measured"
        row["iteration"] = str(iteration)
        rows.append(row)
    return rows


def make_v3_rows(mode="top_k", measured_samples=2):
    rows = make_rows(measured_samples)
    for row in rows:
        row.update(
            {
                "schema_version": "3",
                "mode": mode,
                "max_path_cost_cells": "0",
                "goal_count": "1",
                "completion": "requested_k_reached",
                "per_goal_complete": "true",
                "per_goal_outcomes": "complete",
                "per_goal_limits": "none",
                "per_goal_completions": "requested_k_reached",
                "per_goal_found_paths": "1",
                "max_cost_bounded_paths": "1000",
                "max_path_points": "100000",
                "max_multi_goal_count": "32",
            }
        )
        if mode == "all_within_length":
            row.update(
                {
                    "k": "0",
                    "max_path_cost_cells": "400",
                    "completion": "cost_bound_exhausted",
                    "per_goal_completions": "cost_bound_exhausted",
                }
            )
        elif mode == "multi_goal":
            row.update(
                {
                    "k": "0",
                    "max_path_cost_cells": "400",
                    "goal_count": "2",
                    "completion": "all_goals_complete",
                    "found_paths": "2",
                    "expected_paths": "2",
                    "expanded_nodes": "2",
                    "total_path_points": "4",
                    "per_goal_complete": "true|true",
                    "per_goal_outcomes": "complete|complete",
                    "per_goal_limits": "none|none",
                    "per_goal_completions": (
                        "cost_bound_exhausted|frontier_exhausted"
                    ),
                    "per_goal_found_paths": "1|1",
                }
            )
    return rows


def write_csv(path, rows, fields=FIELDS_V2):
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        writer.writerows(
            {field: row[field] for field in fields}
            for row in rows
        )


class ProfileSummaryContractTest(unittest.TestCase):

    def setUp(self):
        self.temporary_directory = tempfile.TemporaryDirectory()
        self.directory = pathlib.Path(self.temporary_directory.name)

    def tearDown(self):
        self.temporary_directory.cleanup()

    def run_summary(self, paths, expected_samples=2, *extra_arguments):
        return subprocess.run(
            [
                sys.executable,
                str(SUMMARY),
                "--format",
                "csv",
                "--expected-measured-samples",
                str(expected_samples),
                *extra_arguments,
                *(str(path) for path in paths),
            ],
            check=False,
            capture_output=True,
            text=True,
        )

    def test_complete_case_passes(self):
        profile = self.directory / "complete.csv"
        write_csv(profile, make_rows())

        result = self.run_summary([profile])

        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn(",2,PASS,PASS\n", result.stdout)

    def test_schema_v2_is_reported_as_legacy_top_k(self):
        profile = self.directory / "legacy-v2.csv"
        write_csv(profile, make_rows())

        result = self.run_summary([profile])

        self.assertEqual(result.returncode, 0, result.stderr)
        header, record = result.stdout.splitlines()
        self.assertEqual(
            header.split(",")[-3:], ["measured_samples", "verdict", "acceptance"]
        )
        values = dict(zip(header.split(","), record.split(",")))
        self.assertEqual(values["mode"], "top_k")
        self.assertEqual(values["max_path_cost_cells"], "0.000")
        self.assertEqual(values["goal_count"], "1")

    def test_schema_v3_top_k_case_passes(self):
        profile = self.directory / "top-k-v3.csv"
        write_csv(profile, make_v3_rows(), FIELDS_V3)

        result = self.run_summary([profile])

        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn("open_256,top_k,1,0,1,", result.stdout)

    def test_schema_v3_all_within_length_case_passes(self):
        profile = self.directory / "bounded-v3.csv"
        write_csv(profile, make_v3_rows("all_within_length"), FIELDS_V3)

        result = self.run_summary([profile])

        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn("open_256,all_within_length,0,400,1,", result.stdout)
        self.assertTrue(result.stdout.rstrip().endswith(",2,PASS,PASS"))

    def test_complete_bounded_no_path_certificate_passes(self):
        profile = self.directory / "bounded-no-path-v3.csv"
        rows = make_v3_rows("all_within_length")
        for row in rows:
            row["success"] = "false"
            row["outcome"] = "no_path"
            row["found_paths"] = "0"
            row["expected_paths"] = "0"
            row["expanded_nodes"] = "0"
            row["total_path_points"] = "0"
            row["shortest_cost_cells"] = "0.000"
            row["longest_cost_cells"] = "0.000"
            row["per_goal_outcomes"] = "no_path"
            row["per_goal_found_paths"] = "0"
        write_csv(profile, rows, FIELDS_V3)

        result = self.run_summary([profile])

        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertTrue(result.stdout.rstrip().endswith(",2,PASS,PASS"))

    def test_bounded_path_cost_above_bound_cannot_pass(self):
        profile = self.directory / "bounded-cost-overrun.csv"
        rows = make_v3_rows("all_within_length")
        for row in rows:
            row["longest_cost_cells"] = "401.000"
        write_csv(profile, rows, FIELDS_V3)

        result = self.run_summary([profile])

        self.assertEqual(result.returncode, 1, result.stderr)
        self.assertTrue(result.stdout.rstrip().endswith(",FAIL"))

    def test_schema_v3_multi_goal_case_passes(self):
        profile = self.directory / "multi-goal-v3.csv"
        write_csv(profile, make_v3_rows("multi_goal"), FIELDS_V3)

        result = self.run_summary([profile])

        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn("open_256,multi_goal,0,400,2,2,2,", result.stdout)
        self.assertTrue(result.stdout.rstrip().endswith(",2,PASS,PASS"))

    def test_multi_goal_mixed_path_and_no_path_certificates_pass(self):
        profile = self.directory / "mixed-multi-goal-v3.csv"
        rows = make_v3_rows("multi_goal")
        for row in rows:
            row["found_paths"] = "1"
            row["expected_paths"] = "1"
            row["total_path_points"] = "2"
            row["per_goal_outcomes"] = "complete|no_path"
            row["per_goal_found_paths"] = "1|0"
        write_csv(profile, rows, FIELDS_V3)

        result = self.run_summary([profile])

        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertTrue(result.stdout.rstrip().endswith(",2,PASS,PASS"))

    def test_schema_v3_header_must_be_exact(self):
        profile = self.directory / "missing-v3-field.csv"
        write_csv(profile, make_v3_rows(), FIELDS_V3[:-1])

        result = self.run_summary([profile])

        self.assertEqual(result.returncode, 2)
        self.assertIn("schema v2 or v3", result.stderr)

    def test_bounded_mode_rejects_nonzero_k(self):
        profile = self.directory / "bounded-nonzero-k.csv"
        rows = make_v3_rows("all_within_length")
        for row in rows:
            row["k"] = "1"
        write_csv(profile, rows, FIELDS_V3)

        result = self.run_summary([profile])

        self.assertEqual(result.returncode, 2)
        self.assertIn("requires K=0", result.stderr)

    def test_bound_change_between_repeats_is_rejected(self):
        profile = self.directory / "changing-bound.csv"
        rows = make_v3_rows("all_within_length")
        rows[-1]["max_path_cost_cells"] = "401"
        write_csv(profile, rows, FIELDS_V3)

        result = self.run_summary([profile])

        self.assertEqual(result.returncode, 2)
        self.assertIn("measured records", result.stderr)

    def test_multi_goal_array_length_must_match_goal_count(self):
        profile = self.directory / "short-per-goal-array.csv"
        rows = make_v3_rows("multi_goal")
        for row in rows:
            row["per_goal_found_paths"] = "2"
        write_csv(profile, rows, FIELDS_V3)

        result = self.run_summary([profile])

        self.assertEqual(result.returncode, 2)
        self.assertIn("array length does not match goal_count", result.stderr)

    def test_incomplete_multi_goal_certificate_cannot_pass(self):
        profile = self.directory / "incomplete-goal.csv"
        rows = make_v3_rows("multi_goal")
        for row in rows:
            row["completion"] = "incomplete"
            row["per_goal_complete"] = "true|false"
            row["request_satisfied"] = "false"
            row["search_complete"] = "false"
            row["verdict"] = "INCOMPLETE_max_paths"
            row["acceptance"] = "FAIL"
        write_csv(profile, rows, FIELDS_V3)

        result = self.run_summary([profile])

        self.assertEqual(result.returncode, 1, result.stderr)
        self.assertTrue(result.stdout.rstrip().endswith(",2,INCOMPLETE_max_paths,FAIL"))

    def test_multi_goal_requires_at_least_two_goals(self):
        profile = self.directory / "one-goal-multi.csv"
        rows = make_v3_rows("multi_goal")
        for row in rows:
            row["goal_count"] = "1"
            row["per_goal_complete"] = "true"
            row["per_goal_outcomes"] = "complete"
            row["per_goal_limits"] = "none"
            row["per_goal_completions"] = "cost_bound_exhausted"
            row["per_goal_found_paths"] = "2"
        write_csv(profile, rows, FIELDS_V3)

        result = self.run_summary([profile])

        self.assertEqual(result.returncode, 2)
        self.assertIn("goal_count between 2 and 32", result.stderr)

    def test_max_path_points_must_admit_a_complete_path(self):
        profile = self.directory / "invalid-max-path-points.csv"
        rows = make_v3_rows("all_within_length")
        for row in rows:
            row["max_path_points"] = "1"
        write_csv(profile, rows, FIELDS_V3)

        result = self.run_summary([profile])

        self.assertEqual(result.returncode, 2)
        self.assertIn("max_path_points must be at least 2", result.stderr)

    def test_multi_goal_count_cap_must_admit_all_goals(self):
        profile = self.directory / "too-small-multi-goal-cap.csv"
        rows = make_v3_rows("multi_goal")
        for row in rows:
            row["max_multi_goal_count"] = "1"
        write_csv(profile, rows, FIELDS_V3)

        result = self.run_summary([profile])

        self.assertEqual(result.returncode, 2)
        self.assertIn("must be at least goal_count", result.stderr)

    def test_new_resource_limits_must_use_canonical_integers(self):
        profile = self.directory / "noncanonical-new-limit.csv"
        rows = make_v3_rows()
        for row in rows:
            row["max_multi_goal_count"] = "032"
        write_csv(profile, rows, FIELDS_V3)

        result = self.run_summary([profile])

        self.assertEqual(result.returncode, 2)
        self.assertIn("max_multi_goal_count is not canonically encoded", result.stderr)

    def test_resource_limit_change_between_repeats_cannot_pass(self):
        profile = self.directory / "changing-resource-provenance.csv"
        rows = make_v3_rows()
        rows[-1]["max_path_points"] = "99999"
        write_csv(profile, rows, FIELDS_V3)

        result = self.run_summary([profile])

        self.assertEqual(result.returncode, 1, result.stderr)
        self.assertTrue(result.stdout.rstrip().endswith(",2,PASS,FAIL"))

    def test_reported_path_points_cannot_exceed_v3_resource_limit(self):
        profile = self.directory / "path-points-over-v3-limit.csv"
        rows = make_v3_rows("multi_goal")
        for row in rows:
            row["max_path_points"] = "3"
        write_csv(profile, rows, FIELDS_V3)

        result = self.run_summary([profile])

        self.assertEqual(result.returncode, 1, result.stderr)
        self.assertTrue(result.stdout.rstrip().endswith(",FAIL"))

    def test_header_only_input_is_not_silently_ignored(self):
        complete = self.directory / "complete.csv"
        empty = self.directory / "header-only.csv"
        write_csv(complete, make_rows())
        write_csv(empty, [])

        result = self.run_summary([complete, empty])

        self.assertEqual(result.returncode, 2)
        self.assertIn("no profiling records were found", result.stderr)

    def test_truncated_measurements_are_rejected(self):
        profile = self.directory / "truncated.csv"
        write_csv(profile, make_rows(measured_samples=1))

        result = self.run_summary([profile], expected_samples=2)

        self.assertEqual(result.returncode, 2)
        self.assertIn("expected 2 measured records, got 1", result.stderr)

    def test_failure_verdict_cannot_summarize_as_pass(self):
        profile = self.directory / "failed-verdict.csv"
        rows = make_rows()
        for row in rows:
            row["verdict"] = "FAIL_CONTRACT"
        write_csv(profile, rows)

        result = self.run_summary([profile])

        self.assertEqual(result.returncode, 1, result.stderr)
        self.assertIn("FAIL_CONTRACT,FAIL\n", result.stdout)

    def test_path_count_mismatch_cannot_summarize_as_pass(self):
        profile = self.directory / "count-mismatch.csv"
        rows = make_rows()
        for row in rows:
            row["expected_paths"] = "2"
        write_csv(profile, rows)

        result = self.run_summary([profile])

        self.assertEqual(result.returncode, 1, result.stderr)
        self.assertIn("1,2,", result.stdout)
        self.assertTrue(result.stdout.rstrip().endswith(",FAIL"))

    def test_request_satisfied_must_match_terminal_semantics(self):
        profile = self.directory / "request-semantics.csv"
        rows = make_rows()
        for row in rows:
            row["request_satisfied"] = "false"
        write_csv(profile, rows)

        result = self.run_summary([profile])

        self.assertEqual(result.returncode, 1, result.stderr)
        self.assertTrue(result.stdout.rstrip().endswith(",FAIL"))

    def test_expected_count_must_match_the_fixed_scenario(self):
        profile = self.directory / "scenario-count.csv"
        rows = make_rows()
        for row in rows:
            row["found_paths"] = "2"
            row["expected_paths"] = "2"
            row["request_satisfied"] = "false"
        write_csv(profile, rows)

        result = self.run_summary([profile])

        self.assertEqual(result.returncode, 1, result.stderr)
        self.assertTrue(result.stdout.rstrip().endswith(",FAIL"))

    def test_one_case_cannot_be_spliced_across_input_files(self):
        first_half = self.directory / "first-half.csv"
        second_half = self.directory / "second-half.csv"
        rows = make_rows()
        write_csv(first_half, rows[:2])
        write_csv(second_half, rows[2:])

        result = self.run_summary([first_half, second_half])

        self.assertEqual(result.returncode, 2)
        self.assertIn("must come from one input file", result.stderr)

    def test_missing_rss_sample_cannot_summarize_as_pass(self):
        profile = self.directory / "missing-rss.csv"
        rows = make_rows()
        for row in rows:
            row["rss_available"] = "false"
        write_csv(profile, rows)

        result = self.run_summary([profile])

        self.assertEqual(result.returncode, 1, result.stderr)
        self.assertTrue(result.stdout.rstrip().endswith(",FAIL"))

    def test_nonempty_paths_must_contain_at_least_two_points_each(self):
        profile = self.directory / "missing-points.csv"
        rows = make_rows()
        for row in rows:
            row["total_path_points"] = "0"
        write_csv(profile, rows)

        result = self.run_summary([profile])

        self.assertEqual(result.returncode, 1, result.stderr)
        self.assertTrue(result.stdout.rstrip().endswith(",FAIL"))

    def test_malformed_path_digest_is_rejected(self):
        profile = self.directory / "bad-digest.csv"
        rows = make_rows()
        for row in rows:
            row["path_digest"] = "not-a-v2-digest"
        write_csv(profile, rows)

        result = self.run_summary([profile])

        self.assertEqual(result.returncode, 2)
        self.assertIn("invalid path digest", result.stderr)

    def test_wall_time_must_cover_map_and_plan_phases(self):
        profile = self.directory / "impossible-timing.csv"
        rows = make_rows()
        for row in rows:
            row["wall_time_ms"] = "1.0"
        write_csv(profile, rows)

        result = self.run_summary([profile])

        self.assertEqual(result.returncode, 2)
        self.assertIn("wall timing is shorter", result.stderr)

    def test_process_hwm_cannot_decrease_within_a_case(self):
        profile = self.directory / "decreasing-hwm.csv"
        rows = make_rows()
        for row, hwm in zip(rows, ("4096", "1024", "2048")):
            row["process_hwm_kib_after_plan"] = hwm
        write_csv(profile, rows)

        result = self.run_summary([profile])

        self.assertEqual(result.returncode, 2)
        self.assertIn("process HWM decreases", result.stderr)

    def test_path_points_cannot_exceed_the_tree_bound(self):
        profile = self.directory / "excessive-points.csv"
        rows = make_rows()
        for row in rows:
            row["total_path_points"] = "999"
        write_csv(profile, rows)

        result = self.run_summary([profile])

        self.assertEqual(result.returncode, 1, result.stderr)
        self.assertTrue(result.stdout.rstrip().endswith(",FAIL"))

    def test_records_must_remain_in_producer_order(self):
        profile = self.directory / "out-of-order.csv"
        rows = make_rows()
        write_csv(profile, [rows[2], rows[0], rows[1]])

        result = self.run_summary([profile])

        self.assertEqual(result.returncode, 2)
        self.assertIn("producer phase order", result.stderr)

    def test_k_must_use_the_producer_encoding(self):
        profile = self.directory / "noncanonical-k.csv"
        rows = make_rows()
        rows[1]["k"] = "01"
        write_csv(profile, rows)

        result = self.run_summary([profile])

        self.assertEqual(result.returncode, 2)
        self.assertIn("canonically encoded", result.stderr)

    def test_self_crossing_policy_must_match_the_runner(self):
        profile = self.directory / "wrong-policy.csv"
        rows = make_rows()
        for row in rows:
            row["allow_self_crossing"] = "true"
        write_csv(profile, rows)

        result = self.run_summary([profile])

        self.assertEqual(result.returncode, 2)
        self.assertIn("unexpected self-crossing policy", result.stderr)

    def test_standard_matrix_rejects_a_partial_case_set(self):
        profile = self.directory / "one-case.csv"
        write_csv(profile, make_rows())

        result = self.run_summary(
            [profile], 2, "--require-standard-matrix"
        )

        self.assertEqual(result.returncode, 2)
        self.assertIn("standard scenario matrix mismatch", result.stderr)


if __name__ == "__main__":
    unittest.main()
