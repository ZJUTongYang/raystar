#!/usr/bin/env bash

# Run the workspace tests without losing the detailed result report when a
# test fails.  GitHub's shell steps otherwise stop before colcon test-result,
# leaving only a generic exit-code annotation.

set -uo pipefail

test_status=0
colcon test "$@" || test_status=$?

result_status=0
result_report="$(colcon test-result --verbose 2>&1)" || result_status=$?
printf '%s\n' "${result_report}"

if ((test_status != 0 || result_status != 0)); then
  if [[ "${GITHUB_ACTIONS:-}" == "true" ]]; then
    # Workflow command data must escape percent signs and line endings.  Cap
    # the source text so one pathological linter report cannot exceed the
    # annotation service's payload limit.
    annotation="${result_report:0:12000}"
    annotation="${annotation//'%'/'%25'}"
    annotation="${annotation//$'\r'/'%0D'}"
    annotation="${annotation//$'\n'/'%0A'}"
    printf '::error title=colcon test failure::%s\n' "${annotation}"
  fi
  exit 1
fi
