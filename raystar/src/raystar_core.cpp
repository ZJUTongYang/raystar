#include <raystar/raystar_core.h>
#include <chrono>
#include <algorithm>
#include <unordered_set>
#include <cmath>

namespace raystar
{

static std::pair<bool, bool> pnpoly(const std::vector<std::pair<double, double>>& ver,
  double testx, double testy)
{
  bool rightray = false, leftray = false;
  for (size_t i = 0, j = ver.size() - 1; i < ver.size(); j = i++) {
    double cross = (ver[j].first - ver[i].first) * (testy - ver[i].second) /
      (ver[j].second - ver[i].second) + ver[i].first;
    if (((ver[i].second > testy) != (ver[j].second > testy)) && (testx < cross))
      rightray = !rightray;
    if (((ver[i].second > testy) != (ver[j].second > testy)) && (testx > cross))
      leftray = !leftray;
  }
  if (rightray == leftray)
    return {rightray, false};
  return {false, true};
}

Node::Node(const Polymap* /*pMap*/, int Nindex, double start_x, double start_y,
  double Gcost, double Hcost,
  const std::vector<std::pair<double, double>>& visibility_region,
  const std::vector<std::pair<int, int>>& topo_V)
  : Nindex_(Nindex), seed_({start_x, start_y}), start_angle_(0), end_angle_(2 * M_PI),
    parent_index_(-1), start_o_({0.0, 0.0}), end_o_({0.0, 0.0}),
    as_a_child_left_gap_(false), Gcost_(Gcost), Hcost_(Hcost)
{
  seed_ = {start_x, start_y};
  local_shortest_path_.emplace_back(static_cast<int>(start_x), static_cast<int>(start_y));
  path_node_index_.emplace_back(Nindex_);
  V_.assign(visibility_region.begin(), visibility_region.end());
  topo_V_.assign(topo_V.begin(), topo_V.end());
}

Node::Node(const Polymap* /*pMap*/, int Nindex, double seed_x, double seed_y,
  double Gcost, double Hcost, int parent_index,
  const std::vector<std::pair<double, double>>& visibility_region,
  const std::vector<std::pair<int, int>>& topo_V)
  : Nindex_(Nindex), seed_({seed_x, seed_y}), start_angle_(0), end_angle_(0),
    parent_index_(parent_index), start_o_({0.0, 0.0}), end_o_({0.0, 0.0}),
    as_a_child_left_gap_(false), Gcost_(Gcost), Hcost_(Hcost)
{
  seed_ = {seed_x, seed_y};
  V_.assign(visibility_region.begin(), visibility_region.end());
  topo_V_.assign(topo_V.begin(), topo_V.end());
}

void Node::generateChild(const Polymap* pMap)
{
  std::vector<double> theta_list(V_.size());
  for (size_t i = 0; i < V_.size(); ++i)
  {
    theta_list[i] = atan2(V_[i].second - seed_.second, V_[i].first - seed_.first);
    theta_list[i] = start_angle_ + normalize_angle_positive(theta_list[i] - start_angle_);
  }

  std::vector<bool> is_a_gap(V_.size(), false);
  std::vector<bool> is_a_left_gap(V_.size(), false);
  std::vector<int> valid_child_indices;

  std::pair<int, int> topo_V_i, topo_V_next;
  for (size_t i = 0; i < V_.size() - 1; ++i)
  {
    size_t next = i + 1;
    double angle_diff = normalize_angle(theta_list[next] - theta_list[i]);
    constexpr double threshold_squared = 0.0001 * 0.0001;
    if (angle_diff * angle_diff < threshold_squared)
    {
      topo_V_i = topo_V_[i];
      topo_V_next = topo_V_[next];
      if (!pMap->areConsecutive(topo_V_next, topo_V_i))
      {
        is_a_gap[i] = true;
        double disi = (V_[i].first - seed_.first) * (V_[i].first - seed_.first) +
          (V_[i].second - seed_.second) * (V_[i].second - seed_.second);
        double disnext = (V_[next].first - seed_.first) * (V_[next].first - seed_.first) +
          (V_[next].second - seed_.second) * (V_[next].second - seed_.second);
        if (disi > disnext)
          is_a_left_gap[i] = true;
        valid_child_indices.emplace_back(i);
      }
    }
  }

  for (auto i : valid_child_indices)
  {
    size_t next = (i + 1) % V_.size();
    topo_V_i = topo_V_[i];
    topo_V_next = topo_V_[next];
    if (is_a_left_gap[i])
    {
      C_.emplace_back(Child(Nindex_, -1,
        static_cast<int>(V_[next].first), static_cast<int>(V_[next].second), true));

      C_.back().start_angle_ = normalize_angle(theta_list[next]);
      auto next_obs = pMap->getNextObs(topo_V_next);
      double contour_angle = atan2(next_obs.second - V_[next].second,
        next_obs.first - V_[next].first);
      C_.back().end_angle_ = C_.back().start_angle_ +
        normalize_angle_positive(contour_angle - C_.back().start_angle_);

      C_.back().c_obs_index_ = topo_V_next.first;
      C_.back().c_ver_index_ = topo_V_next.second;
      C_.back().o_obs_index_ = topo_V_i.first;
      C_.back().o_ver_index_ = topo_V_i.second;
      C_.back().o_ = V_[i];
    }
    else
    {
      C_.emplace_back(Child(Nindex_, -1,
        static_cast<int>(V_[i].first), static_cast<int>(V_[i].second), false));

      auto prev_obs = pMap->getPrevObs(topo_V_i);
      double contour_angle = atan2(prev_obs.second - V_[i].second,
        prev_obs.first - V_[i].first);
      C_.back().start_angle_ = contour_angle;
      C_.back().end_angle_ = contour_angle +
        normalize_angle_positive(theta_list[i] - contour_angle);

      C_.back().c_obs_index_ = topo_V_i.first;
      C_.back().c_ver_index_ = topo_V_i.second;
      C_.back().o_obs_index_ = topo_V_next.first;
      C_.back().o_ver_index_ = topo_V_next.second;
      C_.back().o_ = V_[next];
    }
  }

  for (size_t i = 0; i < C_.size(); ++i)
    C_[i].Cindex_ = static_cast<int>(i);

  for (auto& child : C_)
  {
    child.c_gcost_ = Gcost_ + std::hypot(seed_.first - child.c_.first,
      seed_.second - child.c_.second);
  }
}

void RaystarCore::outlineMap(std::vector<uint8_t>& costarr, int nx, int ny)
{
  for (int i = 0; i < nx; i++) {
    costarr[i] = 1;
    costarr[(ny - 1) * nx + i] = 1;
  }
  for (int i = 0; i < ny; i++) {
    costarr[i * nx] = 1;
    costarr[i * nx + nx - 1] = 1;
  }
}

void RaystarCore::clearRobotCell(std::vector<uint8_t>& data, int nx, int start_x, int start_y)
{
  data[start_y * nx + start_x] = 0;
}

void RaystarCore::getScopedVisibilityRegion(
  Polymap& theMap, Candidate& the_child,
  std::vector<std::pair<double, double>>& visibility_region,
  std::vector<std::pair<int, int>>& topo_V)
{
  visibility_region.clear();
  topo_V.clear();

  int parent_index = the_child.Nindex_;
  int child_index = the_child.Cindex_;
  auto new_source_point = N_[parent_index].C_[child_index].c_;
  std::vector<std::pair<double, double>> fullV;
  std::vector<std::pair<int, int>> full_topoV;
  theMap.getVisibilityRegion(new_source_point.first, new_source_point.second, fullV, full_topoV);

  std::pair<double, double> start_obs, end_obs;
  std::pair<int, int> start_obs_topo, end_obs_topo;
  double start_angle, end_angle;

  if (N_[parent_index].C_[child_index].is_a_left_gap_)
  {
    start_obs = N_[parent_index].C_[child_index].o_;
    end_obs = theMap.getNextObs({N_[parent_index].C_[child_index].c_obs_index_,
      N_[parent_index].C_[child_index].c_ver_index_});
    start_obs_topo = {N_[parent_index].C_[child_index].o_obs_index_,
      N_[parent_index].C_[child_index].o_ver_index_};
    end_obs_topo = theMap.locateVertex(end_obs.first, end_obs.second);
    start_angle = atan2(start_obs.second - new_source_point.second,
      start_obs.first - new_source_point.first);
    end_angle = atan2(end_obs.second - new_source_point.second,
      end_obs.first - new_source_point.first);
  }
  else
  {
    start_obs = theMap.getPrevObs({N_[parent_index].C_[child_index].c_obs_index_,
      N_[parent_index].C_[child_index].c_ver_index_});
    end_obs = N_[parent_index].C_[child_index].o_;
    start_obs_topo = theMap.locateVertex(start_obs.first, start_obs.second);
    end_obs_topo = {N_[parent_index].C_[child_index].o_obs_index_,
      N_[parent_index].C_[child_index].o_ver_index_};
    start_angle = atan2(start_obs.second - new_source_point.second,
      start_obs.first - new_source_point.first);
    end_angle = atan2(end_obs.second - new_source_point.second,
      end_obs.first - new_source_point.first);
  }
  end_angle = start_angle + normalize_angle_positive(end_angle - start_angle);

  double theta;
  for (size_t i = 0; i < fullV.size(); ++i)
  {
    theta = atan2(fullV[i].second - new_source_point.second,
      fullV[i].first - new_source_point.first);
    theta = start_angle + normalize_angle(theta - start_angle);
    if (theta >= start_angle - 0.0000001 && theta <= end_angle + 0.0000001)
    {
      visibility_region.emplace_back(fullV[i]);
      topo_V.emplace_back(full_topoV[i]);
    }
  }

  int loc = static_cast<int>(std::find_if(visibility_region.begin(), visibility_region.end(),
    [&](const auto& a) { return a.first == start_obs.first && a.second == start_obs.second; })
    - visibility_region.begin());

  if (loc == static_cast<int>(visibility_region.size()))
  {
    visibility_region.insert(visibility_region.begin(), start_obs);
    topo_V.insert(topo_V.begin(), theMap.locateVertex(start_obs.first, start_obs.second));
  }
  else
  {
    visibility_region.erase(visibility_region.begin(), visibility_region.begin() + loc);
    topo_V.erase(topo_V.begin(), topo_V.begin() + loc);
  }

  loc = static_cast<int>(std::find_if(visibility_region.begin(), visibility_region.end(),
    [&](const auto& a) { return a.first == end_obs.first && a.second == end_obs.second; })
    - visibility_region.begin());

  if (loc == static_cast<int>(visibility_region.size()))
  {
    visibility_region.emplace_back(end_obs);
    topo_V.emplace_back(theMap.locateVertex(end_obs.first, end_obs.second));
  }
  else
  {
    visibility_region.erase(visibility_region.begin() + loc + 1, visibility_region.end());
    topo_V.erase(topo_V.begin() + loc + 1, topo_V.end());
  }
}

PlanResult RaystarCore::plan(const GridMap& grid_map,
  int start_x, int start_y, int goal_x, int goal_y,
  int K, bool allow_self_crossing)
{
  PlanResult result;

  GridMap work_map = grid_map;
  clearRobotCell(work_map.data, work_map.width, start_x, start_y);
  outlineMap(work_map.data, work_map.width, work_map.height);

  auto map_start_time = std::chrono::high_resolution_clock::now();
  auto theMap = std::make_shared<Polymap>(work_map, start_x, start_y, goal_x, goal_y);
  auto map_end_time = std::chrono::high_resolution_clock::now();
  result.map_time_ms = std::chrono::duration_cast<std::chrono::microseconds>(
    map_end_time - map_start_time).count() / 1000.0;

  if (!theMap->solution_exist_) {
    result.success = false;
    result.message = "No path exists between start and goal";
    result.polymap = theMap;
    return result;
  }

  auto planner_start_time = std::chrono::high_resolution_clock::now();

  auto comp = [](const Candidate& a, const Candidate& b) { return a.Fcost_ > b.Fcost_; };

  Q_.clear();
  N_.clear();
  Q_.emplace_back(Candidate(-1, -1,
    std::hypot(start_x - goal_x, start_y - goal_y)));

  while (!Q_.empty())
  {
    Candidate best_candidate = Q_[0];
    std::pop_heap(Q_.begin(), Q_.end(), comp);
    Q_.pop_back();

    int parent_index = best_candidate.Nindex_;
    int child_index = best_candidate.Cindex_;

    if (parent_index == -1)
    {
      std::vector<std::pair<double, double>> Vtemp;
      std::vector<std::pair<int, int>> topo_Vtemp;
      theMap->getVisibilityRegion(start_x, start_y, Vtemp, topo_Vtemp);
      N_.emplace_back(Node(theMap.get(), 0, start_x, start_y, 0.0,
        best_candidate.Fcost_, Vtemp, topo_Vtemp));
      N_.back().generateChild(theMap.get());
    }
    else
    {
      auto new_source_point = N_[parent_index].C_[child_index].c_;
      int new_node_index = static_cast<int>(N_.size());

      if (!allow_self_crossing)
      {
        int c_obs = N_[parent_index].C_[child_index].c_obs_index_;
        int c_ver = N_[parent_index].C_[child_index].c_ver_index_;

        std::unordered_set<long long> visited_vertices;
        for (int walk = parent_index; walk > 0; walk = N_[walk].parent_index_)
        {
          int pi = N_[walk].parent_index_;
          if (pi < 0) break;
          for (auto& ch : N_[pi].C_)
          {
            if (ch.c_.first == N_[walk].seed_.first &&
              ch.c_.second == N_[walk].seed_.second)
            {
              visited_vertices.insert(
                static_cast<long long>(ch.c_obs_index_) * 100000LL + ch.c_ver_index_);
              break;
            }
          }
        }

        if (visited_vertices.count(static_cast<long long>(c_obs) * 100000LL + c_ver) > 0)
          continue;
      }

      std::vector<std::pair<double, double>> Vtemp;
      std::vector<std::pair<int, int>> topo_Vtemp;
      getScopedVisibilityRegion(*theMap, best_candidate, Vtemp, topo_Vtemp);

      N_.emplace_back(Node(theMap.get(), new_node_index,
        new_source_point.first, new_source_point.second,
        N_[parent_index].C_[child_index].c_gcost_,
        N_[parent_index].C_[child_index].c_hcost_,
        parent_index, Vtemp, topo_Vtemp));

      N_.back().local_shortest_path_.assign(
        N_[parent_index].local_shortest_path_.begin(),
        N_[parent_index].local_shortest_path_.end());
      N_.back().local_shortest_path_.emplace_back(new_source_point);
      N_.back().path_node_index_.emplace_back(new_node_index);
      N_.back().start_angle_ = N_[parent_index].C_[child_index].start_angle_;
      N_.back().end_angle_ = N_[parent_index].C_[child_index].end_angle_;

      N_.back().generateChild(theMap.get());
    }

    for (auto& child : N_.back().C_)
    {
      child.c_hcost_ = std::hypot(child.c_.first - goal_x, child.c_.second - goal_y);
      Q_.emplace_back(Candidate(child.Nindex_, child.Cindex_, child.c_gcost_ + child.c_hcost_));
      std::push_heap(Q_.begin(), Q_.end(), comp);
    }

    auto b = pnpoly(N_.back().V_, goal_x, goal_y);
    if (b.first || b.second)
    {
      auto locally_shortest_path = N_.back().local_shortest_path_;
      locally_shortest_path.emplace_back(goal_x, goal_y);
      double new_path_length = N_.back().Gcost_ +
        std::hypot(N_.back().seed_.first - goal_x, N_.back().seed_.second - goal_y);

      result.path_solutions.emplace_back(
        PathSolution(locally_shortest_path, new_path_length, N_.back().path_node_index_));
      if (static_cast<int>(result.path_solutions.size()) >= K)
        break;
    }
  }

  auto planner_end_time = std::chrono::high_resolution_clock::now();
  result.plan_time_ms = std::chrono::duration_cast<std::chrono::microseconds>(
    planner_end_time - planner_start_time).count() / 1000.0;

  result.success = !result.path_solutions.empty();
  if (!result.success)
    result.message = "Planning completed but no path found";
  result.polymap = theMap;
  return result;
}

}  // namespace raystar
