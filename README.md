# Ray*: k-Shortest Non-homotopic Path Planning for ROS 2

[![Raystar ROS 2 CI](https://github.com/ZJUTongYang/raystar/actions/workflows/ci.yml/badge.svg?branch=ros2)](https://github.com/ZJUTongYang/raystar/actions/workflows/ci.yml?query=branch%3Aros2)

Ray* computes up to `K` topologically distinct, locally shortest paths on a 2D
occupancy grid. This `ros2` branch provides a standalone ROS 2 planner node, a
cancellable cached-map Action, a compatibility Service, structured completion
diagnostics, bounded visualizations, and an RViz2 Panel.

The planning API itself does not depend on Navigation2. The bundled demo uses
Nav2's map server and lifecycle manager only to publish and activate its test
map.

## Highlights

- Continuous start and goal coordinates, with integer obstacle turning
  vertices retained by the geometric planner.
- Multiple non-homotopic paths ordered by nondecreasing metric cost.
- Cached `nav_msgs/OccupancyGrid` Action workflow, so each goal carries a map
  identity instead of another complete map.
- Cooperative Action cancellation and explicit planning, map, path, debug,
  response, and visualization budgets.
- Machine-readable requested/found/returned counts, completion flags, limits,
  and timing fields.
- Deterministic Core, ROS integration, independent geometry-oracle, RViz
  plugin, sanitizer, coverage, static-analysis, and relocated-install tests.

## Packages

| Package | Purpose |
|---|---|
| `raystar` | Planner Core, ROS 2 node, launch file, bundled map, RViz configuration, profiling runner, and tests |
| `raystar_interfaces` | Cached-map Action, compatibility Service, map identity, path/debug records, and structured result messages |
| `raystar_rviz_plugins` | RViz2 Panel for creating and cancelling cached-map Action requests |

## Supported platforms

| Support level | Ubuntu | ROS 2 | Compiler baseline | Geometry |
|---|---|---|---|---|
| Minimum | 22.04 | Humble | GCC 11, C++17 | CGAL 5.4 with GMP/MPFR |
| Current LTS lane | 24.04 | Jazzy | GCC 13, C++17 | CGAL 5.6.1 with GMP/MPFR |

The RViz Panel requires Qt 5.15. Later CGAL releases may build when they
satisfy the `>= 5.4` source contract, but releases outside the table are not
currently CI-backed compatibility claims. See
[raystar/COMPATIBILITY.md](raystar/COMPATIBILITY.md) for the authoritative
platform, dependency, and relocated-install contract.

## Build from source

The following example uses ROS 2 Humble. Use `jazzy` on Ubuntu 24.04.

```bash
source /opt/ros/humble/setup.bash

mkdir -p ~/raystar_ws/src
cd ~/raystar_ws/src
git clone --branch ros2 --single-branch \
  https://github.com/ZJUTongYang/raystar.git

cd ~/raystar_ws
rosdep update
rosdep install --from-paths src --ignore-src --rosdistro humble -r -y

colcon build --symlink-install \
  --cmake-args -DCMAKE_BUILD_TYPE=RelWithDebInfo
source install/setup.bash
```

If `rosdep` has not been initialized on the machine, run `sudo rosdep init`
once before `rosdep update`. Source both `/opt/ros/<distro>/setup.bash` and the
workspace's `install/setup.bash` in every new terminal.

## Quick start

The bundled launch starts and activates the test map server, starts Raystar,
and opens RViz2:

```bash
ros2 launch raystar raystar_demo.launch.py
```

The supplied RViz configuration already loads the map and MarkerArray
displays and contains one `raystar_rviz_plugins/RaystarPanel`. After the panel
has received the map identity, set the start, goal, and `K`, then press
**Plan**.

Use another map with:

```bash
ros2 launch raystar raystar_demo.launch.py \
  map_yaml:=/absolute/path/to/map.yaml
```

For a headless deployment that already publishes an OccupancyGrid:

```bash
ros2 launch raystar raystar_demo.launch.py \
  start_map_server:=false \
  start_rviz:=false \
  map_topic:=/map
```

The launch file also supports namespaces, a custom Action name, resource-limit
overrides, and optional exposure of the full-map compatibility Service. See
[raystar/LAUNCH.md](raystar/LAUNCH.md) for the complete launch guide.

## ROS interfaces

These are the names resolved from the default node name `raystar`:

| Interface | Default name | Purpose |
|---|---|---|
| Action | `/raystar/plan_paths` | `raystar_interfaces/action/PlanRaystarPaths`; cancellable planning against the cached map identity |
| Service | `/raystar/get_raystar_paths` | `raystar_interfaces/srv/GetRaystarPaths`; compatibility request carrying a complete map |
| Topic | `/raystar/map_status` | `raystar_interfaces/msg/MapStatus`; transient-local identity and metadata for the cached map |
| Topic | `/raystar/non_homotopic_paths` | Complete path visualization snapshot |
| Topic | `/raystar/poly_obstacles` | Simplified polygon-obstacle snapshot |
| Topic | `/raystar/debug_tree` | Bounded search-tree debug snapshot |
| Topic | `/raystar/cdt` | Constrained-triangulation edge snapshot |

The normal client sequence is:

1. Wait for `/raystar/map_status`.
2. Copy its `map_id` into a `PlanRaystarPaths` goal.
3. Send start, goal, `K`, and the request policy flags.
4. Inspect both the Action terminal state and `result_info`.
5. Treat `message` as human-readable diagnostics only.

The bundled launch disables the full-map compatibility Service by default.
Directly starting `raystar_node` retains its compatibility default unless
`enable_legacy_map_service` is overridden. New clients should use the Action,
because a ROS Service request cannot cancel an in-flight plan.

## Input and result contract

Important input rules are:

- The input type is `nav_msgs/OccupancyGrid`, not raw
  `nav2_msgs/Costmap` bytes.
- `map.data` values must be `-1` or in `[0, 100]`. `-1` is the only unknown
  encoding and is traversable only when `allow_unknown=true`.
- Known values at or above `occupied_threshold` are occupied; its default is
  `99`. Lower known values have no path penalty because Raystar is a binary
  geometric planner, not a cost-aware planner.
- The map must have finite values, positive resolution, origin `z=0`, and an
  identity origin rotation. Translated x/y origins are supported, but rotated
  grids are rejected.
- Start and goal must be finite, use exactly the map frame, lie inside the map,
  and be strict interior points of reachable free space. The outermost cell
  ring is reserved as the geometric boundary.
- Raystar does not perform TF transforms.

Returning fewer than `K` paths is not inherently an error. For example,
`STATUS_FEWER_PATHS` with `search_complete=true` means the search frontier was
fully exhausted and fewer admissible topological paths exist.
`STATUS_PARTIAL_SEARCH` instead means a timeout or node limit interrupted the
search, while `STATUS_PARTIAL_OUTPUT` means output budgets omitted one or more
Core paths. Clients should use the structured status, completeness flags,
limits, and requested/found/returned counts rather than parsing text.

Planning capacity is one and requests are not queued. Timeout and cancellation
are cooperative: one CGAL/STL primitive already in progress must return before
the stop request can be observed. The complete parameter, QoS, visualization,
and result contract is in [raystar/README.md](raystar/README.md).

## Multi-scenario acceptance and profiling

An optional deterministic runner exercises the production `RaystarCore`
through six fixed scenarios and `K=1,3,10,50`. The runner is disabled in normal
builds and enabled with `-DRAYSTAR_BUILD_PROFILING=ON`.

Every invocation verifies:

- the input map digest is unchanged;
- completion fields, limits, timings, and expanded-node counts are internally
  consistent;
- continuous endpoints are preserved and all path coordinates/costs are
  finite;
- path cost equals independently reconstructed polyline length;
- paths are unique and sorted by nondecreasing cost;
- independent collision and open-polyline self-intersection oracles accept
  every returned path, including paths returned by an incomplete search;
- every solution has a valid, acyclic root-to-leaf ancestry chain; and
- outcome, limit, found count, path digest, and expanded nodes are deterministic
  across the first invocation, warm-ups, and measured repeats.

### Scenario contract

| Scenario | Grid | Purpose | Expected paths |
|---|---:|---|---:|
| `open_256` | 256×256 | Direct-path and open-map construction baseline | `min(K, 1)` |
| `single_obstacle_256` | 256×256 | Two routes around one central 64×64 obstacle | `min(K, 2)` |
| `narrow_gate_256` | 256×256 | Four-cell gate plus outer detours | `min(K, 5)` |
| `dense_lattice_192` | 192×192 | Search through a deterministic 7×7 obstacle lattice | `K` |
| `large_open_1024` | 1024×1024 | Large-map geometry-construction baseline | `min(K, 1)` |
| `bundled_testmap_50` | 50×50 | Compact irregular-obstacle Core regression scene | `K` |

The bundled PGM scenario uses its file byte-row order directly, matching the
existing Core regression fixture. It does not apply the vertical image
conversion performed by `nav2_map_server`, so it should not be interpreted as
an Action round-trip against the map server.

### Recorded environment and protocol

| Field | Recorded value |
|---|---|
| Date | 2026-07-30 |
| Planner source baseline | `0b61cf5e24f0` (the profiling/documentation change does not modify Core) |
| Build | `RelWithDebInfo`, GCC 11.4.0, C++17 |
| OS / ROS 2 | Ubuntu 22.04.5, Linux 6.8.0-51, ROS 2 Humble |
| Geometry libraries | CGAL 5.4, GMP 6.2.1, MPFR 4.1.0 |
| CPU | Intel Core i9-14900KF, 32 logical CPUs, `performance` governor |
| RAM | 62 GiB |
| Scheduling | No CPU affinity; non-isolated development workstation |
| Per-case process | 1 first invocation, 3 unreported warm-ups, 20 measured repeats |
| Planner policy | `allow_self_crossing=false` |
| Limits | `max_nodes=10000`, `planning_timeout=5000 ms`, `max_k=100` |
| Percentiles | Nearest rank over the 20 measured repeats |

Each `(scenario, K)` used a fresh process. Across 24 cases, all 576 Core
invocations passed, including the 72 unreported warm-ups. The runner emitted
504 first/measured records; all 504 were `PASS`, all returned paths passed the
independent geometry oracles, no resource limit was reached, and every repeat
matched the first invocation's path digest and expanded-node count.

### Recorded results

`Found / expected` uses the scenario contract, not the requested `K`. Thus an
open map returning one path for `K=50` is a successful, complete exhaustion.

| Scenario | K | Found / expected | First wall ms | Wall p50 / p95 ms | Map p50 / p95 ms | Plan p50 / p95 ms | Expanded nodes | Process HWM MiB | Acceptance |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---|
| `open_256` | 1 | 1 / 1 | 181.482 | 175.801 / 212.535 | 175.752 / 212.500 | 0.025 / 0.038 | 1 | 6.484 | **PASS** |
| `open_256` | 3 | 1 / 1 | 180.758 | 176.789 / 187.569 | 176.753 / 187.535 | 0.026 / 0.029 | 1 | 6.234 | **PASS** |
| `open_256` | 10 | 1 / 1 | 172.794 | 177.052 / 190.290 | 177.013 / 190.248 | 0.025 / 0.029 | 1 | 6.484 | **PASS** |
| `open_256` | 50 | 1 / 1 | 172.831 | 177.676 / 198.539 | 177.640 / 198.499 | 0.025 / 0.036 | 1 | 6.484 | **PASS** |
| `single_obstacle_256` | 1 | 1 / 1 | 300.017 | 276.749 / 320.727 | 276.337 / 320.316 | 0.405 / 0.434 | 4 | 6.738 | **PASS** |
| `single_obstacle_256` | 3 | 2 / 2 | 281.104 | 284.213 / 319.339 | 283.640 / 318.769 | 0.554 / 0.670 | 9 | 5.988 | **PASS** |
| `single_obstacle_256` | 10 | 2 / 2 | 287.004 | 286.771 / 305.177 | 286.171 / 304.320 | 0.586 / 0.837 | 9 | 6.738 | **PASS** |
| `single_obstacle_256` | 50 | 2 / 2 | 275.519 | 283.794 / 325.146 | 283.232 / 324.575 | 0.563 / 0.702 | 9 | 6.387 | **PASS** |
| `narrow_gate_256` | 1 | 1 / 1 | 369.641 | 385.252 / 409.473 | 384.954 / 409.172 | 0.276 / 0.367 | 1 | 6.730 | **PASS** |
| `narrow_gate_256` | 3 | 3 / 3 | 396.340 | 384.533 / 421.474 | 383.738 / 420.659 | 0.789 / 1.065 | 7 | 6.234 | **PASS** |
| `narrow_gate_256` | 10 | 5 / 5 | 416.771 | 396.085 / 428.229 | 394.378 / 426.477 | 1.635 / 2.424 | 33 | 6.242 | **PASS** |
| `narrow_gate_256` | 50 | 5 / 5 | 377.427 | 400.555 / 432.175 | 398.602 / 430.507 | 1.648 / 2.565 | 33 | 6.488 | **PASS** |
| `dense_lattice_192` | 1 | 1 / 1 | 1221.824 | 1172.343 / 1249.475 | 1162.060 / 1238.804 | 9.165 / 12.220 | 9 | 6.000 | **PASS** |
| `dense_lattice_192` | 3 | 3 / 3 | 1224.689 | 1203.102 / 1251.517 | 1154.954 / 1207.667 | 41.971 / 54.897 | 67 | 6.500 | **PASS** |
| `dense_lattice_192` | 10 | 10 / 10 | 1242.277 | 1215.812 / 1254.757 | 1153.241 / 1190.483 | 57.401 / 70.949 | 97 | 7.750 | **PASS** |
| `dense_lattice_192` | 50 | 50 / 50 | 1316.722 | 1305.582 / 1378.042 | 1144.558 / 1229.389 | 146.604 / 170.371 | 367 | 10.750 | **PASS** |
| `large_open_1024` | 1 | 1 / 1 | 3007.371 | 3015.035 / 3092.842 | 3014.932 / 3092.691 | 0.028 / 0.041 | 1 | 29.730 | **PASS** |
| `large_open_1024` | 3 | 1 / 1 | 2965.516 | 2981.001 / 3025.851 | 2980.780 / 3025.732 | 0.028 / 0.039 | 1 | 30.527 | **PASS** |
| `large_open_1024` | 10 | 1 / 1 | 3054.373 | 3015.433 / 3062.531 | 3015.313 / 3062.425 | 0.028 / 0.035 | 1 | 30.367 | **PASS** |
| `large_open_1024` | 50 | 1 / 1 | 3080.781 | 3049.434 / 3219.909 | 3049.315 / 3219.789 | 0.029 / 0.047 | 1 | 29.957 | **PASS** |
| `bundled_testmap_50` | 1 | 1 / 1 | 40.643 | 38.217 / 46.753 | 35.354 / 41.998 | 2.824 / 3.025 | 6 | 5.000 | **PASS** |
| `bundled_testmap_50` | 3 | 3 / 3 | 45.951 | 40.799 / 48.479 | 35.459 / 40.952 | 5.245 / 7.159 | 12 | 5.000 | **PASS** |
| `bundled_testmap_50` | 10 | 10 / 10 | 63.734 | 52.097 / 66.994 | 35.097 / 39.942 | 16.889 / 26.657 | 50 | 5.750 | **PASS** |
| `bundled_testmap_50` | 50 | 50 / 50 | 104.836 | 98.394 / 109.998 | 34.990 / 39.254 | 61.096 / 74.161 | 337 | 7.750 | **PASS** |

The open-map rows show that polygon/CDT construction, rather than tree search,
dominates this Core baseline. Increasing `K` mainly affects true multi-route
scenes: for example, dense-lattice plan p50 rises from `9.165 ms` at `K=1` to
`146.604 ms` at `K=50`, while its deterministic tree grows from 9 to 367
nodes. These coarse phases identify where time is spent, but they do not by
themselves justify a particular geometry-index or cache rewrite; that requires
sub-phase counters and before/after semantic regression evidence.

### Reproduce the record

Build the optional tools from the workspace root:

```bash
colcon build --symlink-install \
  --cmake-args \
    -DCMAKE_BUILD_TYPE=RelWithDebInfo \
    -DRAYSTAR_BUILD_PROFILING=ON \
    -DBUILD_TESTING=ON
source install/setup.bash
```

Run each tuple in a separate process so its process HWM is attributable to one
case:

```bash
profile_dir=$(mktemp -d /tmp/raystar-profile.XXXXXX)
profile_bin=$(ros2 pkg prefix raystar)/lib/raystar/raystar_profile

for scenario in \
  open_256 single_obstacle_256 narrow_gate_256 \
  dense_lattice_192 large_open_1024 bundled_testmap_50
do
  for k_value in 1 3 10 50
  do
    "$profile_bin" \
      --scenario "$scenario" \
      --k-values "$k_value" \
      --warmups 3 \
      --iterations 20 \
      > "$profile_dir/${scenario}_K${k_value}.csv"
  done
done

ros2 run raystar raystar_profile_summary \
  --format markdown \
  --process-isolated-memory \
  "$profile_dir"/*_K*.csv
```

Both tools return nonzero when validation, determinism, schema, or acceptance
fails. A quick CTest smoke is also registered whenever profiling is enabled:

```bash
ctest --test-dir build/raystar \
  -R test_profile_smoke --output-on-failure
```

### Measurement semantics and limitations

- `map_time_ms` measures `Polymap::create()`, primarily polygon/CDT geometry
  construction. It is not ROS map-reception time.
- `plan_time_ms` measures endpoint validation and tree search after Polymap
  construction.
- `wall_time_ms` surrounds the complete `RaystarCore::plan()` call. It excludes
  the runner's map digest and independent oracle checks.
- Every call rebuilds its work map and Polymap. “First” and “subsequent” do not
  mean map-cache miss/hit, and the first-call value does not include process
  startup or claim a machine-cold cache.
- The Core runner excludes OccupancyGrid conversion, Action scheduling, DDS
  transport, response serialization, ROS path interpolation, MarkerArray
  construction, and RViz rendering. Measure Action client round-trip time for
  end-to-end deployment latency.
- Process HWM is sampled after planning and is the process-wide high-water
  mark, including the executable, shared libraries, input scenario, Core, and
  earlier runner work in that same tuple. Because HWM is monotonic, measured
  rows can also retain memory peaks from a preceding independent-oracle check
  in the same tuple process. It is not planner-only allocation.
- These values are a reference baseline for the recorded machine, source, and
  build. They are neither hard real-time guarantees nor cross-hardware
  performance claims.

## Tests and CI

Run the normal test suite from the workspace root:

```bash
colcon build --symlink-install \
  --cmake-args -DBUILD_TESTING=ON
colcon test --event-handlers console_cohesion+
colcon test-result --verbose
```

GitHub Actions covers:

- complete Humble/Jammy and Jazzy/Noble builds and tests;
- profiling build and smoke execution;
- AddressSanitizer and UndefinedBehaviorSanitizer;
- native coverage collection;
- clang-format, clang-tidy, CMake lint, and XML contracts; and
- a moved non-symlink install, exported target consumer, pluginlib load,
  installed profiling tools, node startup, and headless launch.

## Documentation

| Document | Contents |
|---|---|
| [raystar/README.md](raystar/README.md) | Detailed planner behavior, parameters, interfaces, result semantics, and RViz behavior |
| [raystar/LAUNCH.md](raystar/LAUNCH.md) | Demo/manual launch, custom maps, namespaces, and deployment settings |
| [raystar/COMPATIBILITY.md](raystar/COMPATIBILITY.md) | Supported platforms, dependency ownership, rosdep keys, and relocation contract |
| [PlanRaystarPaths.action](raystar_interfaces/action/PlanRaystarPaths.action) | Cached-map cancellable Action schema |
| [GetRaystarPaths.srv](raystar_interfaces/srv/GetRaystarPaths.srv) | Full-map compatibility Service schema |
| [PlanningResultInfo.msg](raystar_interfaces/msg/PlanningResultInfo.msg) | Status, completeness, count, limit, and timing fields |

## Citation

If Ray* contributes to your work, please cite:

```bibtex
@inproceedings{Yang2024Tree,
  author    = {Yang, Tong and Huang, Li and Wang, Yue and Xiong, Rong},
  title     = {Tree-based Representation of Locally Shortest Paths for
               2D k-Shortest Non-homotopic Path Planning},
  booktitle = {2024 IEEE International Conference on Robotics and Automation
               (ICRA)},
  year      = {2024}
}
```

## License

Ray* is released under the [MIT License](LICENSE).

Copyright (c) 2024-2026 Tong Yang.
