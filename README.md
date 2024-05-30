# Ray*: Efficient Calculation of 2D Non-homotopic Shortest Paths
Ray* is a C++ ROS implementation of a novel planner to solve the k-shortest non-homotopic path planning (k-SNPP) problem. It is particularly useful in scenarios where a robot needs to compute alternative routes for travel, a group of mobile robots requires multiple non-homotopic routes to avoid conflicts and collisions, or a tethered robot must select suitable topological routes to move without violating maximum tether length constraints. 

The package is compatible with the ```global_planner``` package in ROS Navigation: It visualizes the k-paths as markers and reports the shortest path to the ```move_base``` framework. Designed for grid-based maps, Ray* offers an efficient solution, making it a suitable replacement for the ```global_planner```. 

![Demonstration of Some Tests](https://github.com/ZJUTongYang/raystar/blob/main/doc/raystar_ros_demo.mp4)

## Citation
If you use Ray* in your project, please cite the following paper: 
 ```bibtex
@article{Yang2024Tree, 
author={Tong Yang, Li Huang, Yue Wang, and Rong Xiong}, 
journal={Proceedings of the IEEE International Conference on Robotic and Automation (ICRA) 2024}, 
title={Tree-based Representation of Locally Shortest Paths for 2D k-Shortest Non-homotopic Path Planning},
year={2024}, 
volume={}, 
number={}, 
pages={1-7}, 
}
```

## Setup

### Requirements

(1) CGAL (Must be built from source instead of using ```apt-get```. Version CGAL-5.6.1 is recommended.)

(2) Other ROS-related packages required by ```global_planner``` in ROS Navigation, such as ```move_base```, ```costmap_2d```, etc. 

### Installation
(1) Clone this repository
(2) Build it with ```catkin_make```

### Usage
To use Ray* with ROS Navigation stack, replace the default global planner setting with Ray*: 
(1) Open the ```move_base``` launch file. 
(2) Find the parameter for the global planner and change it from ```global_planner/GlobalPlanner``` to ```raystar/Raystar```. 
(3) Run the ROS Navigation stack as you normally would.

## Interface

### Input:

The number of required non-homotopic locally shortest path, ```K``` (default: ```K = 5```).

### Output

(1) ```non_homotopic_paths```(type: ```visualization_msgs::MarkerArray```): A topic that visualises all the k paths.
(2) ```poly_obstacles```(type: ```visualization_msgs::MarkerArray```): A topic that visualises polygonal obstacles (for debugging purposes).

## License
Ray* is licensed under the MIT license. 

