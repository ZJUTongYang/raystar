# Ray*: 2D Non-homotopic Shortest Path Planner (ROS2)

Ray* is a C++ ROS2 implementation of a planner that solves the k-shortest non-homotopic path planning (k-SNPP) problem. Given a grid map, start/goal positions, and a number K, it computes K topologically distinct paths.

This is the ROS2 port of the [original Ray*](https://github.com/ZJUTongYang/raystar) for ROS1, redesigned as a standalone ROS2 service without Navigation2 dependency.

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

- ROS2 Humble Hawksbill
- CGAL >= 5.6.1 (tested with 5.6.1, compatible with 6.1.1)
- C++17 compiler

## Building

```bash
cd ros2_version
colcon build --symlink-install
source install/setup.bash
```

## Running

```bash
ros2 run raystar raystar_node
```

The node exposes:
- **Service**: `~/get_raystar_paths` — call with map, start, goal, K to get paths
- **Topic**: `~/non_homotopic_paths` (`visualization_msgs/MarkerArray`) — visualization of all K paths
- **Topic**: `~/poly_obstacles` (`visualization_msgs/MarkerArray`) — polygonal obstacles for debugging

## Service Interface

```
raystar_interfaces/srv/GetRaystarPaths

Request:
  nav_msgs/OccupancyGrid    map
  geometry_msgs/PoseStamped start
  geometry_msgs/PoseStamped goal
  int32                     k                    # number of non-homotopic paths
  bool                      allow_self_crossing   # allow paths that wind around obstacles
  bool                      allow_unknown         # treat unknown cells (-1) as free

Response:
  bool                      success               # true if any path found
  string                    message               # error message if failed
  nav_msgs/Path[]           paths                 # resulting paths (world coords)
  float64[]                 path_costs            # cost of each path
```

## RViz2 Visualization

1. Add `MarkerArray` display, topic: `/raystar/non_homotopic_paths`
2. Add `MarkerArray` display, topic: `/raystar/poly_obstacles`
3. Add "Raystar Panel" from the Panels menu for interactive control

## Running Tests

```bash
colcon build --cmake-args -DBUILD_TESTING=ON
colcon test
colcon test-result --verbose
```

## License

MIT
