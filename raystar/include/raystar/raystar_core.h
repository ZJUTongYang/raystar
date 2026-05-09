#pragma once

#include <vector>
#include <utility>
#include <functional>
#include <raystar/coordinate_utils.h>
#include <raystar/polymap.h>

namespace raystar
{

class Candidate
{
public:
  int Nindex_;
  int Cindex_;
  double Fcost_;

  Candidate(int node_index, int child_index, double cost)
    : Nindex_(node_index), Cindex_(child_index), Fcost_(cost) {}

  Candidate(const Candidate&) = default;
};

class Child
{
public:
  int Nindex_;
  int Cindex_;
  double start_angle_;
  double end_angle_;

  std::pair<int, int> c_;
  std::pair<double, double> o_;

  int c_obs_index_;
  int c_ver_index_;
  int o_obs_index_;
  int o_ver_index_;

  bool is_a_left_gap_;

  double c_gcost_;
  double c_hcost_;

  Child(int nindex, int cindex, int cx, int cy, bool is_left_gap)
    : Nindex_(nindex), Cindex_(cindex), start_angle_(0), end_angle_(0),
      c_({cx, cy}), o_({0.0, 0.0}),
      c_obs_index_(-1), c_ver_index_(-1), o_obs_index_(-1), o_ver_index_(-1),
      is_a_left_gap_(is_left_gap), c_gcost_(0), c_hcost_(0) {}
};

class Node
{
public:
  int Nindex_;
  std::pair<int, int> seed_;
  double start_angle_;
  double end_angle_;
  int parent_index_;

  std::pair<double, double> start_o_;
  std::pair<double, double> end_o_;
  bool as_a_child_left_gap_;

  double Gcost_;
  double Hcost_;

  std::vector<Child> C_;

  std::vector<std::pair<double, double>> V_;
  std::vector<std::pair<int, int>> topo_V_;

  std::vector<std::pair<int, int>> local_shortest_path_;
  std::vector<int> path_node_index_;

  Node(const Polymap* pMap, int Nindex, double start_x, double start_y,
    double Gcost, double Hcost,
    const std::vector<std::pair<double, double>>& visibility_region,
    const std::vector<std::pair<int, int>>& topo_V);

  Node(const Polymap* pMap, int Nindex, double seed_x, double seed_y,
    double Gcost, double Hcost, int parent_index,
    const std::vector<std::pair<double, double>>& visibility_region,
    const std::vector<std::pair<int, int>>& topo_V);

  void generateChild(const Polymap* pMap);
};

struct PathSolution
{
  std::vector<std::pair<int, int>> path_;
  double path_cost_;
  std::vector<int> path_node_index_;
  PathSolution(const std::vector<std::pair<int, int>>& the_path,
    double the_path_cost, const std::vector<int>& the_path_node_index)
    : path_cost_(the_path_cost)
  {
    path_.assign(the_path.begin(), the_path.end());
    path_node_index_.assign(the_path_node_index.begin(), the_path_node_index.end());
  }
};

struct PlanResult
{
  bool success = false;
  std::string message;
  std::vector<PathSolution> path_solutions;
  double map_time_ms = 0.0;
  double plan_time_ms = 0.0;
  const Polymap* polymap = nullptr;
};

class RaystarCore
{
public:
  RaystarCore() = default;

  PlanResult plan(const GridMap& grid_map,
    int start_x, int start_y, int goal_x, int goal_y,
    int K, bool allow_self_crossing);

private:
  void outlineMap(std::vector<uint8_t>& costarr, int nx, int ny);
  void clearRobotCell(std::vector<uint8_t>& data, int nx, int start_x, int start_y);

  void getScopedVisibilityRegion(Polymap& theMap,
    Candidate& the_child,
    std::vector<std::pair<double, double>>& visibility_region,
    std::vector<std::pair<int, int>>& topo_V);

  std::vector<Candidate> Q_;
  std::vector<Node> N_;
};

}  // namespace raystar
