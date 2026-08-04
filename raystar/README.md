# Ray*: 2D Non-homotopic Shortest Path Planner (ROS2)

Ray* is a C++ ROS2 planner for topologically distinct locally shortest paths.
Given a grid map and start/goal positions, it computes either the first `K`
paths or all paths whose metric length is at most an inclusive bound. A native
multi-goal mode enumerates bounded paths to several goals while expanding one
shared topology tree. It also implements homotopy-preserving Untethered Path
Shortening (UPS) for directed pairs of tether configurations.

## UPS geometry contract

UPS does not search the graph of triangulation edges. Raystar traces the
reference `alpha_a^-1 * alpha_b` through free CDT faces, cancels immediate
portal reversals, retains repeated face occurrences for winding sleeves, and
runs a funnel over the directed portals. Consequently, the resultant polyline
may cross triangle interiors and turns only where the sleeve requires it.

The C++ API is available as `Polymap::shortenPathWithinHomotopy()`,
`RaystarCore::shortenWithinHomotopy()`, and the explicit-pair batch method
`RaystarCore::shortenConfigurationTransitions()`. The cached-map ROS 2 API is
`/raystar/build_transition_graph`
(`raystar_interfaces/action/BuildRaystarTransitionGraph`). Every successful
transition reports collision-free, homotopy-preserved, and locally-shortest
certificates. A same-configuration pair is a valid zero-length transition.
When a tether configuration originates in a ROS `PathResult`, UPS must receive
its unsampled `topology_path`, not its dense execution/display `path`. The two
have the same endpoints and metric cost, but interpolation roundoff in the
dense path can move a boundary-following segment into an obstacle face.

This is the ROS 2 port of the
[original Ray*](https://github.com/ZJUTongYang/raystar) for ROS 1, redesigned
around standalone ROS 2 Actions with a backward-compatible Service endpoint
and no Navigation2 dependency.

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

The package includes a test map (`test/testmap.pgm` + `test/testmap.yaml`) and
an RViz2 config (`rviz/raystar_test.rviz`). The recommended demonstration is
the one-command launch below; the three-terminal sequence is retained as a
manual bringup and debugging alternative.

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
`goal_set_action_name`,
`transition_action_name`,
`start_map_server`, `start_rviz`, all planner/resource-limit parameters, and
`enable_legacy_map_service`. It defaults to the safer cached-map Action
workflow (`enable_legacy_map_service:=false`). If a namespace or custom Action
name is selected, update the RViz Panel's separate single-goal and multi-goal
Action fields to the corresponding endpoints.

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
- the standard **Publish Point** tool on `/clicked_point` for start/goal capture

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
Invalid requests return `success=false` with an explanatory message. For the
single-goal Action and compatibility Service, a successful response always
contains at least one returned path; a path found by Core but omitted by an
output budget is reported in `message`. The multi-goal aggregate success
contract deliberately permits certified empty results, as documented below.

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

The server applies the following ROS parameters to every request. It reads one
configuration snapshot before map conversion or geometry work:

- `occupied_threshold` (default `99`): known `OccupancyGrid` values at or
  above this value are occupied. It must be between `1` and `100`. This
  occupancy policy is startup-only and read-only while the node runs.
- `max_k` (default `100`): largest accepted value of request `k`.
- `max_cost_bounded_paths` (default `1000`): safety cap on paths retained by an
  `ALL_WITHIN_LENGTH` search. Reaching it before the frontier crosses the
  requested bound reports `LIMIT_MAX_PATHS` and an uncertified partial search.
- `max_multi_goal_count` (default `128`): largest goal array admitted by one
  shared-tree `PlanRaystarGoalSet` request. `max_cost_bounded_paths` applies
  independently to each goal, while `max_path_points` remains request-wide.
- `max_transition_configurations` (default `4096`): largest flattened tether
  configuration array admitted by one `BuildRaystarTransitionGraph` request.
  This is independent of `max_multi_goal_count`.
- `max_transition_pairs` (default `1000`): largest explicit directed pair
  array admitted by one `BuildRaystarTransitionGraph` request. This is
  independent of `max_cost_bounded_paths`.
- `max_nodes` (default `10000`): maximum number of fully constructed search
  nodes, including the root node.
- `planning_timeout_ms` (default `5000`): cooperative deadline for one admitted
  planning request, including map geometry construction, Raystar expansion, or
  UPS portal/funnel processing. It must be positive; the largest `int64`
  millisecond value is reserved by the Core API as its no-deadline sentinel
  and is rejected by the ROS server.
- `max_map_cells` (default `8388608`): maximum `width * height` accepted from
  an `OccupancyGrid`.
- `max_map_bytes` (default `536870912`): maximum map working-set budget.  The
  admission check conservatively budgets 32 bytes per cell for the incoming
  grid, the binary copy, and the known full-grid Core/Polymap scratch arrays;
  this is a conservative admission estimate, not an exact process-RSS limit.
- `max_path_points` (default `100000`): aggregate path-pose budget. Planning
  responses account for both dense and topology paths; UPS admission applies
  it to all input tether poses and UPS output applies it to all returned
  transition poses. A planning path that cannot preserve its complete dense
  and topology representations within the remaining budget is omitted rather
  than partially shortened.
- `max_debug_nodes` (default `0`): maximum search nodes exported in the
  structured `DebugNode[]` response and the debug-tree marker when a request
  explicitly sets `include_debug=true`. A zero budget disables those node
  outputs while retaining planning; a nonzero parameter alone does not opt a
  request into the potentially larger debug payload.
- `max_response_bytes` (default `67108864`): conservative response payload
  budget, including both planning path representations, UPS transitions, and
  debug arrays. The node stops adding paths or debug nodes before this budget
  is exceeded and reports any omission in `message`. The minimum accepted
  value is 1024 bytes; accounting includes CDR alignment headroom and is
  intentionally conservative across RMWs. The same budget is divided across
  the four visualization topics, so polygon, CDT, path, and debug MarkerArrays
  are truncated before their point/marker vectors can grow without bound.

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
  -p max_cost_bounded_paths:=1000 \
  -p max_multi_goal_count:=128 \
  -p max_transition_configurations:=4096 \
  -p max_transition_pairs:=1000 \
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
All 14 integer parameters publish range descriptors. An invalid startup
override makes the node exit before it advertises its Service or Actions
instead of starting with an unusable configuration.

The 12 resource/request-limit parameters from `max_k` through
`max_response_bytes` remain dynamically writable. A pure on-set validator
rejects an invalid update before commit and preserves the previous values. Use
the atomic parameter API when changing several values as one configuration
transaction. The node reads one consistent snapshot before map conversion.
Once that read completes, the request keeps the snapshot; an accepted update
is observed only by requests whose snapshot read happens afterward.
`occupied_threshold` and `path_visualization_republish_period_ms` may be
overridden at startup but reject runtime changes. Keeping the occupancy policy
immutable ensures that GCP and UPS requests handled by one node interpret a
given cached map identity with the same binary occupancy semantics.

A top-K request above `max_k`, a goal set above `max_multi_goal_count`, a UPS
batch above either transition limit, or a request above either map budget is
rejected before its expensive geometry work. If `max_nodes`,
`max_cost_bounded_paths`, or the deadline stops an incomplete Raystar search,
paths found before the limit are returned with an explicit incomplete-search
result; single-goal `success` is true iff at least one path was returned. The
multi-goal Action instead uses aggregate `request_satisfied`, including
certified empty results. UPS timeout is reported separately as
`BuildRaystarTransitionGraph.Result.STATUS_TIMEOUT`; an accepted client
cancellation remains `STATUS_CANCELLED`. If a valid path or debug set exceeds
an output budget, the response contains a bounded subset and reports partial
output plus the omission in `message` (or returns `success=false` when no path
fits). Candidate selection and output-budget admission follow Core's exact
topology-cost ranking; after serialization, each emitted subset is stably
ordered by the final public cost that covers both wire polylines. It is
therefore not promised to be a prefix under public-cost ordering. The queue is
request-local, so a limited request cannot poison the next plan.

The timeout is cooperative rather than a hard real-time preemption boundary:
the planner checks it between geometry and expansion phases, but an individual
CGAL or visibility operation already in progress must return before the
request can stop.

All Actions use the same cooperative mechanism for client cancellation. The
node keeps the executor free while an Action goal runs on one node-owned worker,
so a cancel request can be accepted during a long search. Cancellation is
polled during occupancy conversion, map construction, visibility work, search,
UPS pair processing, and bounded response construction. It still cannot
interrupt a single CGAL or STL primitive that is already executing.

Service requests and all three Action types share one stateful planner and a
capacity-one admission slot. A concurrent Action goal is rejected; a concurrent
Service request returns `success=false` with
`result_info.status=STATUS_BUSY`. Requests are not queued, which prevents an
unbounded backlog of full map copies. The legacy Service does not have a
client-cancel protocol; use an Action when explicit cancel is required.

### Node interface

The node exposes:
- **Action**: `~/plan_paths` — cancellable planning using a `map_id` for the
  immutable `OccupancyGrid` currently cached from the configured `map_topic`
- **Action**: `~/plan_goal_set` — cancellable, exhaustive bounded planning from
  one start to `goals[]`, using one shared Raystar tree
- **Action**: `~/build_transition_graph` — cancellable UPS shortening for an
  explicit ordered array of directed tether-configuration pairs
- **Service**: `~/get_raystar_paths` — compatibility call carrying a complete
  map; it is advertised only when `enable_legacy_map_service=true`
- **Topic**: `~/map_status` (`raystar_interfaces/MapStatus`) — transient-local
  identity, generation, dimensions, startup occupancy threshold, and planning
  semantic versions for the current cached map
- **Topic**: `~/non_homotopic_paths` (`visualization_msgs/MarkerArray`) — visualization of all returned paths (transient-local, with a cached periodic refresh by default)
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
  uint8                     search_mode           # TOP_K or ALL_WITHIN_LENGTH
  int32                     k                     # >0 in TOP_K; 0 in bounded mode
  float64                   max_path_length       # metres; 0 in TOP_K; >0 in bounded mode
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

The two request modes are mutually exclusive. `SEARCH_MODE_TOP_K` requires
`k>0` and `max_path_length=0`; `SEARCH_MODE_ALL_WITHIN_LENGTH` requires `k=0`
and a finite positive `max_path_length`. The length bound is inclusive.

Every returned `PathResult` contains two representations of the same locally
shortest path. `path` is densely sampled for execution, display, and collision
checking. `topology_path` is the unsampled Core polyline containing continuous
endpoints and obstacle turning vertices. Consumers that use the result as a
homotopy reference or stable configuration identity must use `topology_path`.
Both paths are start-to-goal, use the map frame, and have the metric length in
`cost`. More precisely, `cost` is the least binary64 value that
conservatively upper-bounds the mathematical length of both serialized
polylines. Dense interpolation and the grid-to-world transform can make those
two exact lengths differ by a final rounding ULP, so the server finalizes this
public cost after serialization and then stably restores nondecreasing cost
order. Equal public costs retain Core's deterministic topology order.

### Multi-goal Action

`raystar_interfaces/action/PlanRaystarGoalSet` accepts one `start`, ordered
`goals[]`, and an equally sized `max_path_lengths[]` array in metres. The
bounds are positive, finite, and inclusive; they may differ per goal. The
server returns an ordered `GoalPathResult[]`, where `goal_index` and the echoed
pose make duplicate goals unambiguous. Each entry has its own
`PlanningResultInfo`, message, and paths, while `MultiGoalPlanningResultInfo`
reports shared expanded-node and timing data. Debug nodes are request-wide;
each per-goal `result_info` echoes `debug_requested` and the aggregate
`debug_output_complete` certificate.

The Core expands one tree. A candidate is ordered by
`min_i(G + distance(candidate, goal_i) - bound_i)` over active goals. A
strictly positive shared frontier key proves every remaining bound exhausted.
An unreachable goal receives its own `STATUS_NO_PATH` without preventing
reachable goals in the same request from completing. The per-goal
`LIMIT_MAX_PATHS` cap retires only that goal; global node, timeout, cancellation,
path-point, and response-byte limits remain request-wide. To avoid an
`O(frontier * goals)` proof scan after every expansion, a below-cap goal is
retired only by the shared frontier proof (or by final frontier exhaustion).
If a request-wide limit stops another goal first, such a goal is therefore
reported conservatively as partial even when a separate post-hoc scan might
have proved its shorter bound complete; no false completeness is claimed.

The top-level `PlanRaystarGoalSet.Result.success` describes the aggregate
operation, not path presence: it equals aggregate `request_satisfied`. A
complete, fully serialized enumeration therefore succeeds at the ROS Action
transport even when every goal has zero paths. Per-goal `GoalPathResult.success`
means only that `path_results` is non-empty; inspect its `result_info` to
distinguish a certified empty result from partial or failed planning.

The RViz Panel exposes this Action as **Multi-goal: All within lengths**. Its
ordered goal table stores one `(x, y, budget)` row per target and sends one
`PlanRaystarGoalSet` request, so the server expands one shared tree rather than
issuing independent single-goal requests. The per-goal result table keeps
certified empty results visible and reports completeness from
`cost_bound_exhausted && output_complete`. The Action shares cancellation,
cached-map identity, the worker thread, and the capacity-one admission slot
with `~/plan_paths`, the UPS Action, and the legacy Service.

### UPS transition Action

`raystar_interfaces/action/BuildRaystarTransitionGraph` accepts the cached
`map_id`, an ordered `tether_configurations[]` array, and explicit directed
index pairs. Every non-empty configuration must use the exact map frame, run
from one exactly identical base coordinate to its robot endpoint, and preserve
the intended piecewise-linear topology. When the configuration comes from
Raystar planning, pass `PathResult.topology_path`; do not pass the interpolated
`PathResult.path`.

Raystar builds one free-space CDT for the complete batch, traces each reference
`alpha_from^-1 * alpha_to` as a lifted portal sleeve, and emits the funnel
shortening from the source endpoint to the destination endpoint. Pair order and
duplicates are preserved. `triangle_occurrences` is also occurrence-ordered,
so a winding sleeve may repeat the same stable face ID. UPS output is not
restricted to triangulation edges and may cross triangle interiors.

A capacity-one cache retains only a completed, immutable Polymap. Reuse
requires the exact cached-map generation, occupancy policy, tether base, and
canonical endpoint set to match a completed GCP or UPS construction. This
enables the GCP-to-UPS pipeline without weakening cache identity. Any map
update or cache-key mismatch rebuilds the environment; mutable search-tree
state is never cached.

Every admitted planning result also returns a versioned `environment_id`.
This identity is stable across timestamp-only republishes and binds the map
frame and geometry, raw occupancy payload, startup-only `occupied_threshold`,
the request's `allow_unknown` value, and Raystar's occupancy, geometry, and
topology semantic versions. `map_id` continues to identify the exact cached
message snapshot; the two identities serve different purposes.

`BuildRaystarTransitionGraph.Goal.expected_environment_id` is an optional
fail-closed guard. A zero UUID preserves the source behavior of older UPS
clients and derives the current environment from `map_id` and
`allow_unknown`. A non-zero UUID must match exactly before any transition
geometry is built. The UPS result echoes the actual `environment_id` once its
cached map and startup configuration have been resolved. A GCP-to-UPS pipeline
should pass the GCP result identity as this expected value.

`requested_transition_count` is the requested pair count;
`completed_transition_count` is the number of ordered pair records evaluated
and appended, not the number certified. `success=true` and `STATUS_COMPLETE`
require every pair to be collision-free, homotopy-preserving, and locally
shortest. A geometric failure still returns its `HomotopyTransitionResult`,
but the aggregate result is `STATUS_FAILED`. Invalid input, cooperative
deadline expiration, and accepted client cancellation are distinguished by
`STATUS_INVALID_REQUEST`, `STATUS_TIMEOUT`, and `STATUS_CANCELLED`; only the
last uses the ROS Action terminal state `CANCELED`.

`max_transition_configurations` and `max_transition_pairs` are independent UPS
admission limits. They deliberately do not reuse `max_multi_goal_count` or
`max_cost_bounded_paths`. `max_path_points`, `max_response_bytes`, map budgets,
and `planning_timeout_ms` also apply to the UPS request.

UPS Action feedback reports validation, transition-environment preparation,
and pair-shortening stages. `requested_transition_count` stays constant and
`completed_transition_count` is monotonic. Large batches are evenly throttled:
the first and last pair updates are retained and a completed request publishes
at most 100 feedback messages. Feedback is a volatile, non-authoritative live
stream rather than a complete event log: a client may miss an ordered prefix
before its goal handle is registered, or intermediate samples due to QoS
history. It is emitted only while the goal remains active; feedback allocation
or middleware failures do not change planning. The Action Result is the
authoritative terminal progress record.

Before sending a goal, subscribe to `~/map_status` (transient-local), copy the
current `map_id`, and use that identity in the goal. If a newer map arrives
before execution, the server rejects the stale identity instead of silently
planning against a different map. `computeMapId()` is a deterministic
identity hash for routing/version checks; it is not an authentication token.
`computeEnvironmentId()` is likewise a deterministic identity guard, not a
cryptographic proof.

The added fields change the generated ROS interface type hashes. Source-level
UPS compatibility is retained by the zero expected-identity default, but all
Raystar interface producers and consumers in an overlay must be rebuilt
together; do not mix binaries generated from the old and new Actions.

For single- and multi-goal planning, `allow_self_crossing=false` treats the
result as an open polyline: consecutive segments may meet at their shared
construction waypoint, while every intersection between non-consecutive
closed segments is rejected. This includes proper crossings, endpoint or
tangential contact, positive-length collinear overlap, and revisiting an
earlier waypoint. `allow_self_crossing=true` disables that pruning only; it
does not require the planner to manufacture a self-intersecting path or
guarantee that more paths will be returned. Obstacle collision, visibility,
map, and resource-limit checks remain active in either mode.

## Compatibility Service Interface

```
raystar_interfaces/srv/GetRaystarPaths

Request:
  nav_msgs/OccupancyGrid    map                  # values -1 or 0..100; >= occupied_threshold is occupied
  geometry_msgs/PoseStamped start                # finite position in the exact map frame
  geometry_msgs/PoseStamped goal                 # finite position in the exact map frame
  uint8                     search_mode          # top-K or exhaustive within length
  int32                     k                    # >0 in top-K; 0 in bounded mode
  float64                   max_path_length      # metres; inclusive bounded-mode threshold
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

- `requested_path_count` is the requested K (zero in bounded mode). In bounded
  mode, `found_path_count` counts paths certified eligible under the original
  metric bound from their serialized-world topology geometry, before output
  resource filtering; `returned_path_count` equals `path_results.size()`;
- `request_satisfied` means all K paths were returned in top-K mode; in bounded
  mode it means `cost_bound_exhausted && output_complete`;
- `search_complete` means Core found K, exhausted its search frontier, or
  exhausted the requested cost bound;
- `cost_bound_exhausted` is the bounded-search completeness certificate: no
  unexpanded branch can produce a path whose length is at most the bound;
- `output_complete` means every authoritative metric-eligible path was
  included in the ROS response. Candidates admitted only by the padded Core
  search superset and then certified outside the original metric bound are not
  omissions;
- `limits_reached` is a bitmask for timeout, max-nodes, max-paths, path-point,
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
search. If `found_path_count > returned_path_count`, ROS output limits or a
publication-certificate failure omitted one or more metric-eligible paths.
`message` is for people and must not be parsed by clients.

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

The Panel has separate single-goal and multi-goal Action fields, defaulting to
`/raystar/plan_paths` and `/raystar/plan_goal_set`. They accept relative or
absolute ROS names, so a saved RViz configuration can target namespaced
deployments. Editing either endpoint cancels any in-flight request before its
client is rebuilt. The endpoints, map and clicked-point topics, start/single
goal, search mode, K, maximum length, ordered multi-goal rows with independent
budgets, and request policy flags are persisted with the RViz configuration.

The three Panel modes are **Single: Top K**, **Single: All within length**, and
**Multi-goal: All within lengths**. The last mode always uses bounded exhaustive
enumeration and disables K. In multi-goal mode, changing **Set all budgets**
immediately replaces every existing row's budget. Individual `Budget (m)` cells
remain editable afterward, so goals can still use different limits. The
displayed per-goal `Complete` column means
`cost_bound_exhausted && output_complete`; a row may therefore be complete
with status `No path` and zero path payloads. The path table is grouped by goal
and shows the expected flattened `path_N` Marker namespace for visual
cross-reference. Partial results explicitly show that markers were not
published. Even for a complete result, Marker output remains
visualization-only and is not a completeness certificate.

For point-and-click input, keep the Panel's clicked-point topic aligned with
RViz's standard **Publish Point** tool (the bundled configuration uses
`/clicked_point`). **Capture start** and the single-goal **Capture goal** are
one-shot. In multi-goal mode, **Capture goals** remains active and appends one
row per click using the current bulk-budget value; toggle it off when finished
and edit individual budgets directly in the table. Clicked points
must already use the cached map frame because the Panel does not perform TF
transforms.

A map/topic change, receipt of a newer map, the 60-second panel timeout, Action
endpoint change, or panel destruction invalidates the displayed request and
asks the server to cancel an accepted goal. If invalidation occurs before the
goal-acceptance response, a late accepted goal is canceled immediately.

## Deterministic profiling runner

The optional `raystar_profile` executable exercises the Core API without a ROS
graph and treats correctness certificates as acceptance criteria. Build it in
Release mode with:

```bash
colcon build --packages-select raystar --cmake-args \
  -DRAYSTAR_BUILD_PROFILING=ON -DCMAKE_BUILD_TYPE=Release
source install/setup.bash
```

The default invocation remains backward-compatible top-K profiling. The other
two modes perform exhaustive inclusive cost-bounded enumeration for one goal
or expand one shared tree for several goals:

```bash
# Top-K (the --mode option may be omitted)
ros2 run raystar raystar_profile --scenario single_obstacle_256 \
  --mode top-k --k-values 1,3,10 --warmups 3 --iterations 20

# All non-homotopic paths to one goal whose Core cost is at most 300
ros2 run raystar raystar_profile --scenario single_obstacle_256 \
  --mode all-within-length --max-path-length 300 \
  --warmups 3 --iterations 20

# The same bound for four goals, using one shared Raystar tree
ros2 run raystar raystar_profile --scenario dense_lattice_192 \
  --mode multi-goal --max-path-length 300 --goal-count 4 \
  --warmups 3 --iterations 20
```

`--max-path-length` is expressed in Core grid-coordinate units (cells), not
metres. Multiply a reported cell cost by the scenario map resolution to obtain
metres. Multi-goal generation is deterministic: goal zero is the catalogued
scenario goal, then free interior cells are scanned from the largest `y` and
largest `x`; their cell centres are selected while the start and duplicate
goals are skipped. `--goal-count` is limited to 2 through 32.

The runner exposes `--max-nodes`, `--timeout-ms`,
`--max-cost-bounded-paths`, `--max-path-points`, and
`--max-multi-goal-count`. A zero timeout intentionally requests an immediate
cooperative timeout. Invalid option combinations exit with status 2; a valid
run that is incomplete, non-deterministic, or fails an independent path,
terminal-state, or scenario-contract check exits with status 1. Thus a partial
result caused by `max_paths`, node, path-point, or timeout limits is never
reported as an accepted exhaustive enumeration.

CSV schema v3 preserves all 39 schema-v2 columns as an exact prefix and appends
`mode`, `max_path_cost_cells`, `goal_count`, `completion`, the five
pipe-separated per-goal certificate arrays, `max_cost_bounded_paths`,
`max_path_points`, and `max_multi_goal_count`. Single-goal rows contain one
item in every per-goal array. A multi-goal aggregate `completion` is
`all_goals_complete` only when every goal independently reports
`cost_bound_exhausted` or `frontier_exhausted`; the per-goal arrays remain the
authoritative evidence.

The installed `raystar_profile_summary` accepts exact schema-v2 files as
legacy top-K input and exact schema-v3 files for all three modes. Its standard
matrix check remains the six catalogued top-K scenarios at K=1,3,10,50:

```bash
ros2 run raystar raystar_profile_summary profile.csv \
  --expected-measured-samples 20 --format markdown
```

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
