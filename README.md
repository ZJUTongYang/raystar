# Ray*: Non-homotopic Path Planning for ROS 2

[![Raystar ROS 2 CI](https://github.com/ZJUTongYang/raystar/actions/workflows/ci.yml/badge.svg?branch=ros2)](https://github.com/ZJUTongYang/raystar/actions/workflows/ci.yml?query=branch%3Aros2)

Ray* finds topologically distinct, locally shortest paths on a 2D occupancy
grid. It supports top-K planning, exhaustive planning within a length budget,
shared-tree multi-goal planning, and homotopy-preserving Untethered Path
Shortening (UPS). This repository provides the planner, ROS 2 interfaces, and
an RViz2 Panel as one versioned project.

## Quick demo with RaystarPanel

Once the workspace has been [built](#build-from-source), source it and launch
the bundled demo from the workspace root. The commands below use Humble;
Jazzy users should source `/opt/ros/jazzy/setup.bash` instead:

```bash
source /opt/ros/humble/setup.bash
source install/setup.bash
ros2 launch raystar raystar_demo.launch.py
```

The launch file starts and activates the bundled map server, starts Raystar,
and opens an RViz configuration containing `RaystarPanel`.

1. Wait until the Panel reports `Map received: ...` and **Plan** becomes
   enabled.
2. Click **Plan**. No values need to be entered for the first run.
3. Inspect the paths in the **Non-Homotopic Paths** display and their costs and
   completion state in the Panel.

The preset request demonstrates **Multi-goal: All within lengths** with start
`(2.3, 3.3)` and two goals:

| Goal | Inclusive length budget |
|---|---:|
| `(6.7, 7.7)` | 7 m |
| `(6.7, 3.3)` | 6 m |

### Try your own request

Choose one of the three Panel modes:

| Panel mode | Request |
|---|---|
| **Single: Top K** | Up to `K` topologically distinct paths to one goal |
| **Single: All within length** | All admissible paths to one goal whose total cost is at most the inclusive **Max length (m)** |
| **Multi-goal: All within lengths** | Bounded paths to several goals, with one inclusive budget per goal and one shared Raystar tree |

Coordinates can be entered directly or selected from the map:

1. Click **Capture start** or **Capture goal**/**Capture goals**.
2. Select RViz's **Publish Point** tool and click the map.
3. Set `K`, **Max length (m)**, or each table-row **Budget (m)** as required.
4. Click **Plan** again.

In multi-goal mode, **Capture goals** continuously appends rows; click it again
to stop. Select the preset rows and click **Remove selected** first when you
want to replace rather than extend the demo goal set. Changing **Set all
budgets (m)** updates every current row, after which individual budgets can be
edited.

If a custom RViz configuration does not contain the Panel, open **Panels → Add
New Panel** and select `raystar_rviz_plugins/RaystarPanel`.

## Features

- Top-K non-homotopic locally shortest paths.
- All admissible topologically distinct, locally shortest paths within an
  inclusive metric length budget.
- Multi-goal bounded planning with independent goal budgets and one shared
  topology tree.
- Homotopy-preserving UPS transition planning over a reusable raw-contour
  constrained Delaunay triangulation. Inputs must be taut by default, while an
  explicit policy admits complete collision-free references that still need
  shortening; result paths may cross triangle interiors.
- Cached-map, cancellable ROS 2 Actions with structured completion and resource
  diagnostics.
- Point-and-click single- and multi-goal workflows in RaystarPanel.

## One repository, three ROS 2 packages

This is one Raystar repository and release, split into three packages for ROS
interface generation and RViz plugin build separation. The packages are built,
versioned, and deployed together; they are not independent products.

| Package | Purpose |
|---|---|
| `raystar` | Planner Core, ROS 2 node, launch files, bundled map, visualization, profiling tools, and tests |
| `raystar_interfaces` | Planning and UPS Actions, compatibility Service, and result messages |
| `raystar_rviz_plugins` | `RaystarPanel` for interactive single- and multi-goal planning |

The planning API itself does not depend on Navigation2. The demo uses Nav2's
map server and lifecycle manager only to publish and activate its test map.

## Supported platforms

| Support level | Ubuntu | ROS 2 | Compiler | CGAL |
|---|---|---|---|---|
| Minimum | 22.04 | Humble | GCC 11, C++17 | 5.4 with GMP/MPFR |
| Current LTS lane | 24.04 | Jazzy | GCC 13, C++17 | 5.6.1 with GMP/MPFR |

`CGAL >= 5.4` is accepted by the build. The two versions in the table are the
CI-verified distribution baselines; newer releases, including CGAL 6.x, are
not yet a CI-backed compatibility claim. See
[COMPATIBILITY.md](raystar/COMPATIBILITY.md) for details.

## Build from source

The following example uses ROS 2 Humble. Replace `humble` with `jazzy` on
Ubuntu 24.04.

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

Run `sudo rosdep init` once if rosdep has not been initialized. In each new
terminal, source both the ROS distribution and workspace overlay before using
Raystar.

### Upgrading from 0.2.x

Raystar 0.3.0 changes ROSIDL wire types. `BuildRaystarTransitionGraph.Goal`
gains `reference_path_policy` and renames `tether_configurations` to
`rooted_reference_paths`, `ConfigurationTransitionPair` becomes
`ReferenceTransitionPair` with `from_reference`/`to_reference` indices, and
the node parameter `max_transition_configurations` becomes
`max_transition_references`. Stop every running Raystar node, RViz
instance, and downstream client before upgrading. Rebuild
`raystar_interfaces`, `raystar`, `raystar_rviz_plugins`, and every downstream
package together in a clean overlay, then restart every process from a fresh
terminal that sources only the 0.3.0 overlay. A new workspace with empty
`build`, `install`, and `log` directories is the safest upgrade path; mixed
0.2.x/0.3.0 build or runtime overlays are not supported.

Version 0.3.0 also advances `geometry_semantics_version` to 2. Environment
identities therefore change even for a byte-identical occupancy grid, and
clients must refresh them from `MapStatus` after restart. UPS now retains raw
reachable contours for reference shortening. Those contours consume the
unused `max_map_bytes` budget, so a map that passes normal planning admission
can still require a larger byte budget for transition construction. See
[COMPATIBILITY.md](raystar/COMPATIBILITY.md) and
[LAUNCH.md](raystar/LAUNCH.md) for the complete migration and resource
contracts.

## Maps and deployment

To run the Panel demo on another Nav2 map YAML:

```bash
ros2 launch raystar raystar_demo.launch.py \
  map_yaml:=/absolute/path/to/map.yaml
```

To use an existing `nav_msgs/OccupancyGrid` publisher without RViz:

```bash
ros2 launch raystar raystar_demo.launch.py \
  start_map_server:=false \
  start_rviz:=false \
  map_topic:=/map
```

Namespaces, endpoint remapping, map requirements, resource parameters, and
manual bringup are documented in [LAUNCH.md](raystar/LAUNCH.md).

## ROS API at a glance

The default node name is `raystar`:

| Interface | Default endpoint | Purpose |
|---|---|---|
| Action | `/raystar/plan_paths` | Cancellable single-goal top-K or all-within-length planning against the cached map |
| Action | `/raystar/plan_goal_set` | Cancellable multi-goal all-within-length planning using one shared tree |
| Action | `/raystar/build_transition_graph` | UPS for explicit directed transitions between rooted reference paths |
| Service | `/raystar/get_raystar_paths` | Full-map compatibility API for legacy clients; disabled by the demo launch |
| Topic | `/raystar/map_status` | Cached-map identity, semantic versions, and environment identities for both `allow_unknown` policies |

An Action client first waits for `/raystar/map_status` and copies its `map_id`.
For result validation, select `environment_id_disallow_unknown` or
`environment_id_allow_unknown` to match the request. A transition client may
also copy that value into `expected_environment_id`.

The complete schemas and their status/completion semantics are documented in
the interface definitions:

- [PlanRaystarPaths.action](raystar_interfaces/action/PlanRaystarPaths.action)
- [PlanRaystarGoalSet.action](raystar_interfaces/action/PlanRaystarGoalSet.action)
- [BuildRaystarTransitionGraph.action](raystar_interfaces/action/BuildRaystarTransitionGraph.action)
- [GetRaystarPaths.srv](raystar_interfaces/srv/GetRaystarPaths.srv)
- [MapStatus.msg](raystar_interfaces/msg/MapStatus.msg)
- [PathResult.msg](raystar_interfaces/msg/PathResult.msg)
- [PlanningResultInfo.msg](raystar_interfaces/msg/PlanningResultInfo.msg)

The most important client-side rules are:

- Input maps are `nav_msgs/OccupancyGrid` messages. They must be axis-aligned;
  Raystar does not apply TF transforms.
- Start and goal frame IDs must exactly match the map frame.
- Length budgets are finite, positive, inclusive, and expressed in metres.
- A resource limit can produce a partial result. Use the structured status and
  completion fields instead of assuming every returned list is exhaustive.
- `PathResult.path` is the execution/display path and is normally densely
  sampled. When a Raystar result becomes a UPS rooted reference, pass its
  unsampled `PathResult.topology_path`.
- `BuildRaystarTransitionGraph.reference_path_policy` defaults to
  `REFERENCE_PATHS_MUST_BE_TAUT`. Use
  `REFERENCE_PATHS_MAY_BE_UNTAUT` only when a complete reference is
  collision-free in the target map but still needs local shortening, such as
  a configuration-space path shortened against a less restrictive obstacle
  map. This policy does not relax frame, endpoint, common-base, or collision
  checks.

## Tests and profiling

Run the normal suite from the ROS 2 workspace root:

```bash
colcon build --symlink-install --cmake-args -DBUILD_TESTING=ON
colcon test --event-handlers console_cohesion+
colcon test-result --verbose
```

CI runs the full Humble/Jammy and Jazzy/Noble suites together with sanitizer,
coverage, lint, plugin-loading, and relocated-install checks.

The optional deterministic profiling runner is disabled by default:

```bash
colcon build --packages-select raystar --cmake-args \
  -DRAYSTAR_BUILD_PROFILING=ON -DCMAKE_BUILD_TYPE=Release
source install/setup.bash
ros2 run raystar raystar_profile --help
```

It supports top-K, single-goal all-within-length, and shared-tree multi-goal
scenarios. Its command-line length values use Core grid units rather than ROS
API metres; consult `--help` before comparing results.

## Documentation

| Document | Contents |
|---|---|
| [Launch guide](raystar/LAUNCH.md) | Panel bringup, custom maps, namespaces, parameters, and deployment |
| [Compatibility contract](raystar/COMPATIBILITY.md) | Supported platforms, CGAL dependencies, rosdep, and install relocation |
| [ROS interfaces](raystar_interfaces/action/PlanRaystarPaths.action) | Authoritative request, result, status, and completion schemas |

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
