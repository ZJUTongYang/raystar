# Ray*: 2D Non-homotopic Shortest Path Planner (ROS2)

Ray* is a C++ ROS2 implementation of a planner that solves the k-shortest non-homotopic path planning (k-SNPP) problem. Given a grid map, start/goal positions, and a number K, it computes K topologically distinct paths.

This is the ROS2 port of the [original Ray*](https://github.com/ZJUTongYang/raystar) for ROS1, redesigned as a standalone ROS2 Action with a backward-compatible Service endpoint and no Navigation2 dependency.

## Citation

```bibtex
@article{Yang2024Tree,
author={Yang, Tong and Huang, Li and Wang, Yue and Xiong, Rong},
journal={Proceedings of the IEEE International Conference on Robotic and Automation (ICRA) 2024},
title={Tree-based Representation of Locally Shortest Paths for 2D k-Shortest Non-homotopic Path Planning},
year={2024},
}
```

## Requirements

- Ubuntu 22.04 with ROS 2 Humble (minimum supported platform), or Ubuntu 24.04
  with ROS 2 Jazzy (current LTS compatibility lane)
- CGAL >= 5.4 with GMP and MPFR support; CI covers CGAL 5.4 and 5.6.1
- Qt 5.15 or newer within the Qt 5 series for the RViz Panel
- C++17 compiler

See [COMPATIBILITY.md](COMPATIBILITY.md) for the authoritative platform matrix,
dependency ownership, rosdep keys, and relocated-install contract. CGAL 6.x is
accepted by the CMake lower-bound check but is not currently a CI-backed
compatibility claim; Qt 6 is not supported.

## Building

```bash
cd ros2_version
colcon build --symlink-install
source install/setup.bash
```

## Running

### Quick Start with the bundled test map

The package includes a test map (`test/testmap.pgm` + `test/testmap.yaml`) and an RViz2 config (`rviz/raystar_test.rviz`). The full launch flow requires **3 terminals**.

**Prerequisites**: install `nav2_map_server` (one-time):

```bash
sudo apt install ros-humble-nav2-map-server
```

For a one-command demo (including lifecycle activation of the bundled map
server and optional RViz2), run:

```bash
ros2 launch raystar raystar_demo.launch.py
```

The launch file accepts `map_yaml`, `namespace`, `map_topic`, `action_name`,
`start_map_server`, `start_rviz`, all planner/resource-limit parameters, and
`enable_legacy_map_service`. It defaults to the safer cached-map Action
workflow (`enable_legacy_map_service:=false`). If a namespace or custom Action
name is selected, update the RViz panel's **Planner Action / Name** field to
that endpoint.

Each terminal must have the ROS2 + workspace environment sourced (add to `~/.bashrc` for convenience):

```bash
source /opt/ros/humble/setup.bash
source <workspace>/install/setup.bash
```

#### Terminal 1: Map Server

```bash
ros2 run nav2_map_server map_server --ros-args \
  -p yaml_filename:=$(ros2 pkg prefix raystar)/share/raystar/test/testmap.yaml \
  -p use_sim_time:=false
```

Then activate the map server (separate terminal or one-liner):

```bash
ros2 lifecycle set /map_server configure
ros2 lifecycle set /map_server activate
```

Verify the map is published:

```bash
ros2 topic echo /map --once --field info.width   # should print 50
```

#### Terminal 2: Ray* Node

```bash
ros2 run raystar raystar_node
```

#### Terminal 3: RViz2

```bash
rviz2 -d $(ros2 pkg prefix raystar)/share/raystar/rviz/raystar_test.rviz
```

The RViz2 config pre-loads:
- **Map** display on `/map`
- **MarkerArray** on `/raystar/non_homotopic_paths` (per-path namespace toggle)
- **MarkerArray** on `/raystar/poly_obstacles`
- one **RaystarPanel** (`raystar_rviz_plugins/RaystarPanel`) for interactive planning

### Using a custom map

Replace the `yaml_filename` in Terminal 1 with your own map YAML. Ray* currently
accepts only axis-aligned 2D occupancy grids: the YAML `origin` yaw must be
`0.0`, and the published `OccupancyGrid` origin must have `z = 0` and an
identity quaternion. Translated map origins are supported.

### Coordinate and frame contract

The Action and compatibility Service validate their ROS coordinate inputs before invoking the grid
algorithm:

- `map.header.frame_id` must be non-empty.
- `start.header.frame_id` and `goal.header.frame_id` must exactly equal the map
  frame and each frame ID is limited to 256 UTF-8 bytes. The node does not
  perform TF transforms.
- Map resolution must be finite and greater than zero. Map origin coordinates
  and start/goal coordinates must also be finite.
- The map origin and start/goal positions must have `z = 0`. Rotated or tilted
  occupancy grids are rejected.
- Start and goal x/y positions must lie inside the map and be strict interiors
  of the reachable free space; points on an obstacle or map boundary are
  rejected. The planner reserves the outermost ring of grid cells as its map
  boundary, so an endpoint whose floored cell lies in that ring is also
  rejected. Their pose orientations are ignored by this 2D planner.

Custom frame names are supported when the map, start, and goal all use the same
name. Returned `Path` messages and visualization markers use that map frame.
Invalid requests return `success=false` with an explanatory message. A
successful response always contains at least one returned path; a path found
by Core but omitted by an output budget is reported in `message`.

Internally, the planner keeps start/goal in continuous grid coordinates while
obstacle turning points remain integer grid vertices.  Returned paths are
sampled along those segments without rounding; therefore intermediate
`PoseStamped` samples may be fractional even though the turning points are
integer vertices.

### OccupancyGrid value contract

`map.data` must contain `-1` or a value in `[0, 100]`.  The value `-1` is the
only unknown encoding and is treated as free only when the request sets
`allow_unknown=true`; values below `-1` and above `100` reject the complete
request.  Known cells are converted to the binary map used by Raystar with the
server's `occupied_threshold`: values greater than or equal to the threshold
are occupied and lower values are free.

The default threshold is `99`, matching a Nav2 costmap after it has been
translated to `nav_msgs/OccupancyGrid`: Nav2 raw costs `253` and `254` become
`99` and `100`, while raw `255` becomes `-1`.  Set the threshold to `1` to
recover the previous conservative policy in which every nonzero known value is
blocked, or to `100` to block only occupancy value `100`.  This conversion
does not make Raystar cost-aware: values below the threshold carry no path
penalty after conversion.  Raw `nav2_msgs/Costmap` messages with `uint8`
values `0..255` are not accepted by this interface.

### Planning configuration and resource limits

The server applies the following ROS parameters to every request.  It reads
one configuration snapshot before converting the map or invoking Core:

- `occupied_threshold` (default `99`): known `OccupancyGrid` values at or
  above this value are occupied. It must be between `1` and `100`.
- `max_k` (default `100`): largest accepted value of request `k`.
- `max_nodes` (default `10000`): maximum number of fully constructed search
  nodes, including the root node.
- `planning_timeout_ms` (default `5000`): cooperative deadline for one
  `RaystarCore::plan()` call, including map geometry construction and tree
  expansion. It must be positive; the largest `int64` millisecond value is
  reserved by the Core API as its no-deadline sentinel and is rejected by the
  ROS server.
- `max_map_cells` (default `8388608`): maximum `width * height` accepted from
  an `OccupancyGrid`.
- `max_map_bytes` (default `536870912`): maximum map working-set budget.  The
  admission check conservatively budgets 32 bytes per cell for the incoming
  grid, the binary copy, and the known full-grid Core/Polymap scratch arrays;
  this is a conservative admission estimate, not an exact process-RSS limit.
- `max_path_points` (default `100000`): maximum total sampled `PoseStamped`
  points across all paths in one response and its cached visualization.  A
  path that cannot preserve all of its integer turning points and interpolation
  samples within the remaining budget is omitted rather than partially
  shortcut.
- `max_debug_nodes` (default `0`): maximum search nodes exported in the
  structured `DebugNode[]` response and the debug-tree marker when a request
  explicitly sets `include_debug=true`. A zero budget disables those node
  outputs while retaining planning; a nonzero parameter alone does not opt a
  request into the potentially larger debug payload.
- `max_response_bytes` (default `67108864`): conservative response payload
  budget, including path poses and debug arrays.  The node stops adding paths
  or debug nodes before this budget is exceeded and reports any omission in
  `message`.  The minimum accepted value is 1024 bytes; accounting includes
  CDR alignment headroom and is intentionally conservative across RMWs.  The
  same budget is divided across the four visualization topics, so polygon,
  CDT, path, and debug MarkerArrays are truncated before their point/marker
  vectors can grow without bound.

The visualization-only parameter
`path_visualization_republish_period_ms` is read when the node starts. Its
default is `2000`; set it to `0` to disable periodic path refresh. RViz Humble
deletes a path namespace when its checkbox is cleared but does not request a
transient-local replay when the namespace is enabled again, so the default
periodically republishes the already-built path MarkerArray. It does not rerun
planning or rebuild path geometry. This parameter is read-only after startup;
restart the node to change its refresh policy.

For example:

```bash
ros2 run raystar raystar_node --ros-args \
  -p occupied_threshold:=99 \
  -p max_k:=20 \
  -p max_nodes:=5000 \
  -p planning_timeout_ms:=3000 \
  -p max_map_cells:=1000000 \
  -p max_map_bytes:=33554432 \
  -p max_path_points:=50000 \
  -p max_debug_nodes:=500 \
  -p max_response_bytes:=33554432 \
  -p path_visualization_republish_period_ms:=2000
```

`occupied_threshold` must be between 1 and 100. All resource limits except
`max_debug_nodes` must be positive; `max_path_points` must be at least 2 and
`max_map_bytes` at least 32 bytes, and `max_response_bytes` at least 1024.
Every parameter publishes an integer range descriptor. An invalid startup
override makes the node exit before it advertises its Service or Action instead
of starting with an unusable configuration.

The nine planning parameters remain dynamically writable. A pure on-set
validator rejects an invalid update before commit and preserves the previous
values. Use the atomic parameter API when changing several values as one
configuration transaction. The node reads one consistent snapshot before map
conversion. Once that read completes, the request keeps the snapshot; an
accepted update is observed only by requests whose snapshot read happens
afterward. The startup-only
`path_visualization_republish_period_ms` parameter rejects runtime changes.

A request above `max_k` or either map budget is rejected before map
construction. If
`max_nodes` or the deadline stops an incomplete search, paths found before the
limit are returned with `success=true` and an explicit incomplete-search
message; if no path was found, `success=false` is returned. If a valid path or
debug set exceeds an output budget, the response contains the bounded prefix
and explains the omission in `message` (or returns `success=false` when no
path fits). The queue is request-local, so a limited request cannot poison the
next plan.

The timeout is cooperative rather than a hard real-time preemption boundary:
the planner checks it between geometry and expansion phases, but an individual
CGAL or visibility operation already in progress must return before the
request can stop.

The Action uses the same cooperative mechanism for client cancellation.  The
node keeps the executor free while an Action goal runs on one node-owned worker,
so a cancel request can be accepted during a long search.  Cancellation is
polled during occupancy conversion, map construction, visibility work, search,
and bounded response construction.  It still cannot interrupt a single CGAL or
STL primitive that is already executing.

Service requests and Action goals share one stateful planner and a capacity-one
admission slot.  A concurrent Action goal is rejected; a concurrent Service
request returns `success=false` with `result_info.status=STATUS_BUSY`. Requests are not queued,
which prevents an unbounded backlog of full map copies.  The legacy Service
does not have a client-cancel protocol; use the Action when explicit cancel is
required.

### Node interface

The node exposes:
- **Action**: `~/plan_paths` — cancellable planning using a `map_id` for the
  immutable `OccupancyGrid` currently cached from the configured `map_topic`
- **Service**: `~/get_raystar_paths` — compatibility call carrying a complete
  map; it is advertised only when `enable_legacy_map_service=true`
- **Topic**: `~/map_status` (`raystar_interfaces/MapStatus`) — transient-local
  identity, generation, dimensions, and frame for the current cached map
- **Topic**: `~/non_homotopic_paths` (`visualization_msgs/MarkerArray`) — visualization of all K paths (transient-local, with a cached periodic refresh by default)
- **Topic**: `~/poly_obstacles` (`visualization_msgs/MarkerArray`) — polygonal obstacles for debugging (transient-local one-shot snapshot)
- **Topic**: `~/debug_tree` (`visualization_msgs/MarkerArray`) — expanded nodes, visibility regions, and costs
- **Topic**: `~/cdt` (`visualization_msgs/MarkerArray`) — constrained and unconstrained CDT edges

Each MarkerArray publication is a complete replacement snapshot: it clears the
previous markers before adding the current state. Every admitted request clears
the previous visualization when execution starts; if that request later fails
or is canceled, the visualization remains cleared. A request rejected as busy
does not disturb the current snapshot. The transient-local publishers retain
the latest snapshot for late RViz subscribers.
Only the path snapshot is periodically republished, and the timer reuses its
bounded immutable MarkerArray instead of repeating interpolation or Marker
construction. This lets a stock RViz MarkerArray display restore a previously
disabled `path_N` namespace after it is checked again without repeatedly
rebuilding the obstacle, CDT, or debug-tree visualizations.

## Action Interface

```
raystar_interfaces/action/PlanRaystarPaths

Goal:
  unique_identifier_msgs/UUID map_id             # identity from ~/map_status
  geometry_msgs/PoseStamped start
  geometry_msgs/PoseStamped goal
  int32                     k
  bool                      allow_self_crossing
  bool                      allow_unknown         # treat only unknown cells (-1) as free
  bool                      include_debug         # opt in to bounded debug output

Result:
  bool                                  success
  raystar_interfaces/PlanningResultInfo result_info
  string                                message
  raystar_interfaces/PathResult[]       path_results
  raystar_interfaces/DebugNode[]        debug_nodes
```

The Action terminal state is authoritative: a client-accepted cancel completes
with Action state `STATUS_CANCELED`, an empty `path_results` array, and
`result_info.status=STATUS_CANCELLED`. An invalid/no-path goal is aborted; a
complete or partial result containing at least one path succeeds. The Action
currently has no periodic feedback payload.

Before sending a goal, subscribe to `~/map_status` (transient-local), copy the
current `map_id`, and use that identity in the goal. If a newer map arrives
before execution, the server rejects the stale identity instead of silently
planning against a different map. `computeMapId()` is a deterministic
identity hash for routing/version checks; it is not an authentication token.

For both APIs, `allow_self_crossing=false` treats the result as an open
polyline: consecutive segments may meet at their shared construction
waypoint, while every intersection between non-consecutive closed segments is
rejected. This includes proper crossings, endpoint or tangential contact,
positive-length collinear overlap, and revisiting an earlier waypoint.
`allow_self_crossing=true` disables that pruning only; it does not require the
planner to manufacture a self-intersecting path or guarantee that more paths
will be returned. Obstacle collision, visibility, map, and resource-limit
checks remain active in either mode.

## Compatibility Service Interface

```
raystar_interfaces/srv/GetRaystarPaths

Request:
  nav_msgs/OccupancyGrid    map                  # values -1 or 0..100; >= occupied_threshold is occupied
  geometry_msgs/PoseStamped start                # finite position in the exact map frame
  geometry_msgs/PoseStamped goal                 # finite position in the exact map frame
  int32                     k                    # number of non-homotopic paths
  bool                      allow_self_crossing   # apply the open-polyline policy documented above
  bool                      allow_unknown         # treat only unknown cells (-1) as free
  bool                      include_debug         # opt in to bounded debug output

Response:
  bool                                  success       # true iff path_results is non-empty
  raystar_interfaces/PlanningResultInfo result_info   # authoritative structured outcome
  string                                message       # bounded human-readable diagnostic only
  raystar_interfaces/PathResult[]       path_results  # each record binds a Path to its metric cost
  raystar_interfaces/DebugNode[]        debug_nodes   # bounded structured search-tree prefix
```

Every returned `Path` pose uses the map frame.  Pose orientation is the identity
quaternion because this planner computes a 2D geometric path; orientation has
no planning meaning.  The Service fields intentionally mirror the Action goal
and result, but Service clients cannot explicitly cancel an in-flight call.

`PlanningResultInfo` is the machine-readable completion contract:

- `requested_path_count` is the requested K, `found_path_count` is the number
  found by Core, and `returned_path_count` equals `path_results.size()`;
- `request_satisfied` is true only when all K requested paths were returned;
- `search_complete` means Core either found K or exhausted its search frontier;
- `output_complete` means every Core path was included in the ROS response;
- `limits_reached` is a bitmask for timeout, max-nodes, path-point,
  response-byte, and cancellation limits;
- `expanded_nodes`, `map_time_ms`, and `plan_time_ms` provide structured
  diagnostics; `debug_output_complete` reports whether all expanded nodes were
  included in `debug_nodes`.

The high-level status distinguishes `STATUS_COMPLETE`, `STATUS_FEWER_PATHS`,
`STATUS_NO_PATH`, `STATUS_PARTIAL_SEARCH`, `STATUS_PARTIAL_OUTPUT`, invalid
input/configuration, busy, canceled, and internal failure. For example,
`returned_path_count < requested_path_count` with `search_complete=true` and
`STATUS_FEWER_PATHS` means fewer paths exist under the current constraints;
the same count with `STATUS_PARTIAL_SEARCH` means a limit interrupted the
search. If `found_path_count > returned_path_count`, ROS output limits omitted
one or more paths. `message` is for people and must not be parsed by clients.

This revision changes the Service and Action ROS type hashes. All server and
client packages must be rebuilt and deployed together; the term "compatibility
Service" refers to retaining a non-cancellable Service workflow, not binary
wire compatibility with clients compiled against the previous IDL.

## RViz2 Visualization

1. Add `MarkerArray` display, topic: `/raystar/non_homotopic_paths`
2. Add `MarkerArray` display, topic: `/raystar/poly_obstacles`
3. Add `raystar_rviz_plugins/RaystarPanel` from the Panels menu for interactive control

The stable plugin lookup ID is `raystar_rviz_plugins/RaystarPanel`. RViz files
saved by an earlier development version with the implicit C++ type ID
`raystar_rviz_plugins::RaystarPanel` must replace that entry with the stable ID.

The panel's **Planner Action / Name** field defaults to
`/raystar/plan_paths`. It accepts relative or absolute ROS names, so a saved
RViz configuration can target namespaced deployments such as
`/robot1/raystar/plan_paths`. Editing the endpoint cancels any request sent to
the old server before rebuilding the client. The Action name, map topic,
start/goal, K, `allow_self_crossing`, and `allow_unknown` values are persisted
with the RViz configuration.

A map/topic change, receipt of a newer map, the 60-second panel timeout, Action
endpoint change, or panel destruction invalidates the displayed request and
asks the server to cancel an accepted goal. If invalidation occurs before the
goal-acceptance response, a late accepted goal is canceled immediately.

## Running Tests

```bash
colcon build --cmake-args -DBUILD_TESTING=ON
colcon test
colcon test-result --verbose
```

The RViz plugin package includes an offscreen pluginlib smoke test that loads,
constructs, and destroys the canonical Panel ID and verifies configuration
round-tripping. The `raystar` package also checks that the bundled RViz config
contains exactly one canonical Panel entry.

The planner regression suite includes a fixed-seed corpus of 24 small
`OccupancyGrid` maps spanning reachable/unreachable cases and both unknown-cell
policies. An independent four-neighbor BFS supplies the reachability oracle;
successful ROS responses are additionally checked with independent collision
and self-intersection oracles, metric-cost reconstruction, and a repeated-call
determinism check. A failure prints the case seed and complete raw grid so the
same counterexample can be reproduced directly.

## License

MIT
