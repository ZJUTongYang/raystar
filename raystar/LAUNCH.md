# Ray* Launch Guide

## Prerequisites (one-time)

```bash
sudo apt install ros-humble-nav2-map-server
```

 Add to `~/.bashrc` (also one-time):

```bash
source /opt/ros/humble/setup.bash
```

Then source the built overlay from the ROS 2 workspace root in every terminal
before launching (one-time per terminal):

```bash
# run this from the workspace root, for example ~/raystar_ws
source install/setup.bash
```

## Launch

### One-command demo

The package also installs a parameterized bringup that starts the bundled map
server, lifecycle manager, Raystar node, and RViz2 together:

```bash
ros2 launch raystar raystar_demo.launch.py
```

Useful overrides include `map_yaml:=/path/to/map.yaml`,
`namespace:=robot1`, `map_topic:=/robot1/map`,
`action_name:=/robot1/raystar/plan_paths`,
`transition_action_name:=/robot1/raystar/build_transition_graph`,
`goal_set_action_name:=/robot1/raystar/plan_goal_set`, `start_rviz:=false`, and
the planner/resource parameters (`max_k`, `max_cost_bounded_paths`,
`max_multi_goal_count`, `max_transition_references`,
`max_transition_pairs`, `max_nodes`, `planning_timeout_ms`, `max_debug_nodes`,
and so on). The default launch intentionally disables the legacy full-map
Service; set `enable_legacy_map_service:=true` only for old clients. When using
a namespace or custom endpoint, set the RViz Panel's persisted single-goal
and multi-goal Action fields to the corresponding `action_name` and
`goal_set_action_name` values.

### Manual three-terminal alternative

#### Terminal 1: Map Server

```bash
ros2 run nav2_map_server map_server --ros-args \
  -p yaml_filename:=$(ros2 pkg prefix raystar)/share/raystar/test/testmap.yaml \
  -p use_sim_time:=false
```

Activate (in any terminal):

```bash
ros2 lifecycle set /map_server configure
ros2 lifecycle set /map_server activate
```

Verify:

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

## Using a custom map

Replace `yaml_filename` in Terminal 1 with your map YAML path. Ray* accepts only
axis-aligned 2D occupancy grids: resolution must be finite and positive, origin
`z` and yaw must be `0.0`, and the origin quaternion must represent the identity
rotation. Translated x/y origins are supported.

The published map must have a non-empty `header.frame_id` no longer than 256
UTF-8 bytes. Planning start and goal poses must use exactly that frame name and
the same 256-byte bound; Ray* does not perform TF transforms. Non-finite
coordinates and out-of-map start/goal positions are rejected with an Action or Service
error instead of being passed to the planner.

`OccupancyGrid.data` must contain `-1` or values in `[0,100]`. The server
defaults to `occupied_threshold=99`: known values `99` and `100` are occupied,
lower known values are free, and `-1` is controlled only by the request's
`allow_unknown` flag. Values outside that range reject the request. Set the
threshold to `1` for the previous conservative all-nonzero policy. This input
is the translated `nav_msgs/OccupancyGrid` representation; raw
`nav2_msgs/Costmap` values `0..255` are not accepted directly.

The server also defaults to the following configuration and resource limits:
`occupied_threshold=99`, `max_k=100`, `max_cost_bounded_paths=1000`,
`max_multi_goal_count=128`,
`max_transition_references=4096`, `max_transition_pairs=1000`,
`max_nodes=10000`,
`planning_timeout_ms=5000`,
`max_map_cells=8388608`, `max_map_bytes=536870912`,
`max_path_points=100000`, `max_debug_nodes=0` (debug output is opt-in), and
`max_response_bytes=67108864`. The visualization refresh parameter
`path_visualization_republish_period_ms` defaults to `2000`. Override them when
launching if a deployment needs different resource bounds or refresh policy:

```bash
ros2 run raystar raystar_node --ros-args \
  -p occupied_threshold:=99 \
  -p max_k:=20 -p max_cost_bounded_paths:=1000 \
  -p max_multi_goal_count:=128 \
  -p max_transition_references:=4096 \
  -p max_transition_pairs:=1000 \
  -p max_nodes:=5000 -p planning_timeout_ms:=3000 \
  -p max_map_cells:=1000000 -p max_map_bytes:=67108864 \
  -p max_path_points:=50000 -p max_debug_nodes:=500 \
  -p max_response_bytes:=33554432 \
  -p path_visualization_republish_period_ms:=2000
```

`occupied_threshold` must be between 1 and 100. All resource limits except
`max_debug_nodes` must be positive; `max_path_points` must be at least 2,
`max_map_bytes` at least 32 bytes, and `max_response_bytes` at least 1024. A
map must satisfy both `max_map_cells` and the fixed `max_map_bytes` admission
estimate of 32 bytes per cell; this check happens before either Core or Polymap
copies the binary grid.

Ordinary planning retains that simplified-contour admission rule. UPS
transition construction additionally charges its unsimplified reachable
contour topology against the unused byte budget before constructing the CDT.
After fixed map admission, its conservative raw-contour ceiling is:

```text
floor((max_map_bytes - 32 * cell_count) / 4096) vertices
```

The 4096-byte charge is a fail-closed complexity estimate for the coexisting
CGAL and standard-library structures, not an exact heap-size claim. A map can
therefore pass normal planning admission but fail UPS reference-shortening
admission. That failure returns `STATUS_INVALID_REQUEST` with a diagnostic
naming the raw contour and `max_map_bytes`; it does not return a partial
transition graph. Increase `max_map_bytes` when a high-perimeter obstacle map
needs a larger transition environment. A goal set above
`max_multi_goal_count`, or a UPS batch above
`max_transition_references`/`max_transition_pairs`, is rejected before its
expensive geometry work. The two UPS limits are independent of the multi-goal
and per-goal path limits. Reaching the node or time limit stops work cleanly;
Raystar planning reports partial-search fields, while UPS reports
`STATUS_TIMEOUT`. Path representations, UPS inputs/outputs, debug nodes, and
the response itself are bounded as well. The same response budget is split
among the four MarkerArray topics to cap polygon/CDT/path/debug point and marker
vectors. A single CGAL/visibility operation cannot be preempted midway, so the
timeout is cooperative rather than a hard real-time deadline.

All 14 integer parameters expose range descriptors. Invalid launch/YAML
overrides make the process exit before the planning endpoints are advertised.
The 12 resource/request-limit parameters from `max_k` through
`max_response_bytes` may be changed at runtime; invalid updates are rejected
without replacing the previous configuration. Use an atomic parameter update
when several values must change together. Each request keeps the one
configuration snapshot it read before map conversion. An accepted update is
observed only by requests whose snapshot read happens afterward.
`occupied_threshold` is startup-only and read-only while the node runs, so GCP
and UPS interpret a given cached map identity with one binary occupancy policy.

`path_visualization_republish_period_ms` is read once at node startup. A
positive value periodically republishes the cached path MarkerArray so that a
stock RViz Humble MarkerArray display can restore an individually re-enabled
`path_N` namespace. This is a visualization refresh, not a new plan or Marker
rebuild. Set it to `0` to rely only on transient-local delivery and eliminate
periodic path traffic. Obstacle, CDT, and debug-tree snapshots are never
periodically rebuilt or republished. The parameter is read-only while the node
is running; restart the node to change it.

New clients should use the cancellable
`/raystar/plan_paths` (`raystar_interfaces/action/PlanRaystarPaths`) Action.
Multi-goal bounded clients use `/raystar/plan_goal_set`
(`raystar_interfaces/action/PlanRaystarGoalSet`); it accepts ordered `goals[]`
and equally sized inclusive `max_path_lengths[]` arrays and expands one shared
tree. UPS clients use `/raystar/build_transition_graph`
(`raystar_interfaces/action/BuildRaystarTransitionGraph`) with an ordered
configuration array and explicit directed index pairs. `reference_path_policy`
defaults to `REFERENCE_PATHS_MUST_BE_TAUT`; this preserves the original strict
requirement for recompiled zero-initialized clients. Set
`REFERENCE_PATHS_MAY_BE_UNTAUT` when each complete reference is collision-free
in the cached target map but must first be shortened there. Raystar validates
and shortens every complete reference before removing common prefixes and
processing pairs. The opt-in mode does not relax the frame, endpoint,
common-base, or collision requirements.

Transition construction uses a separately cached raw-contour environment;
normal planning geometry is simplified and is not reused for this purpose.
This distinction is required when a reference produced in a stricter map, such
as an inflated configuration space, is shortened in a less restrictive
obstacle map. When a reference comes from planning, pass
`PathResult.topology_path`, not the densely sampled `PathResult.path`.
Cancellation reaches the server's cooperative stop token while the executor
remains available to process the cancel request.  The existing
`/raystar/get_raystar_paths` Service remains available for compatibility, but
the ROS Service protocol cannot cancel an in-flight request.  Action and
Service calls share a capacity-one planner; concurrent work is rejected rather
than queued.  Cancellation also remains cooperative: one already-running CGAL
primitive must return before the goal can reach `STATUS_CANCELED`.

The bundled RViz configuration contains exactly one
`raystar_rviz_plugins/RaystarPanel`. It persists separate single-goal and
multi-goal endpoints; set them to names such as
`/robot1/raystar/plan_paths` and `/robot1/raystar/plan_goal_set` for
multi-robot deployments. Changing either field cancels the in-flight goal.
The bundled Panel opens with two valid goals and independent 7 m and 6 m
budgets in **Multi-goal: All within lengths**, so the default shared-tree
demonstration only requires pressing **Plan** after the map arrives. Its
standard RViz **Publish Point** tool and Panel capture buttons can replace the
default start/goals interactively on
`/clicked_point`.

In multi-goal mode, changing **Set all budgets** updates every existing target
row and invalidates the old result before the next **Plan**. Edit individual
`Budget (m)` cells afterward when goals need different limits.

Raystar 0.4.0 renames the transition-graph fields to planner-neutral names
(`tether_configurations` to `rooted_reference_paths`,
`ConfigurationTransitionPair` to `ReferenceTransitionPair`,
`max_transition_configurations` to `max_transition_references`) and changes
the corresponding DDS type hashes. Environment identities and geometry
semantics are unchanged. Apply the same lockstep rebuild-and-restart
procedure as below.

Raystar 0.3.0 changes the `BuildRaystarTransitionGraph` goal and `MapStatus`
ROS types and therefore changes their DDS type hashes. It also advances
`geometry_semantics_version` to 2, so every environment identity changes even
when the occupancy-grid bytes do not. Stop all running nodes, RViz instances,
and clients; perform a clean lockstep rebuild of all three Raystar packages and
downstream consumers; then restart every process from a fresh terminal that
sources only the 0.3.0 overlay. Refresh identities from the new `MapStatus`.
Endpoint compatibility does not provide wire compatibility with 0.2.x
binaries.

## Stopping

`Ctrl+C` in each terminal in reverse order (rviz2 → raystar_node → map_server).
