#pragma once

#include <iostream>
#include <vector>
#include <algorithm>
#include <unordered_map>
#include <raystar/coordinate_utils.h>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Triangulation_face_base_with_info_2.h>

namespace raystar
{

namespace constrained_delaunay_triangulation
{
  typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
  typedef CGAL::Constrained_Delaunay_triangulation_2<K, CGAL::Default, CGAL::No_constraint_intersection_tag> CDT;
  typedef CDT::Point Point;

  struct BungiuEdge
  {
    std::pair<double, double> prev_pos;
    std::pair<double, double> next_pos;
    std::pair<int, int> topo_prev;
    std::pair<int, int> topo_next;
    double limit_prev;
    double limit_next;
    std::pair<int, int> limit_prev_pos;
    std::pair<int, int> limit_next_pos;
    bool is_bd;
  };
}

class Obs
{
public:
  std::vector<std::pair<int, int>> detailed_ordered_vertices_;
};

class Polymap
{
public:
  int xsize_;
  int ysize_;
  std::vector<uint8_t> data_;

  std::vector<Obs> obs_;

  bool solution_exist_;

  Polymap(const GridMap& grid_map, int start_x, int start_y, int goal_x, int goal_y);
  ~Polymap() = default;

  bool getPolyObstacles(int start_x, int start_y, int goal_x, int goal_y);

  void getVisibilityRegion(int start_x, int start_y,
    std::vector<std::pair<double, double>>& visibility_region,
    std::vector<std::pair<int, int>>& topo_V);

  void simplifyPolyObstacles(int start_x, int start_y, int goal_x, int goal_y);

  std::pair<int, int> locateVertex(int x, int y) const;

  inline bool areConsecutive(const std::pair<int, int>& prev, const std::pair<int, int>& next) const
  {
    return (prev.first == next.first) &&
      (next.second - prev.second + static_cast<int>(obs_[prev.first].detailed_ordered_vertices_.size())) % static_cast<int>(obs_[prev.first].detailed_ordered_vertices_.size()) == 1;
  }

  inline std::pair<int, int> getPrevObs(const std::pair<int, int>& curr) const
  {
    return obs_[curr.first].detailed_ordered_vertices_[(curr.second - 1 + obs_[curr.first].detailed_ordered_vertices_.size()) % obs_[curr.first].detailed_ordered_vertices_.size()];
  }

  inline std::pair<int, int> getNextObs(const std::pair<int, int>& curr) const
  {
    return obs_[curr.first].detailed_ordered_vertices_[(curr.second + 1) % obs_[curr.first].detailed_ordered_vertices_.size()];
  }

  bool calculateVisibilityRegion(int x, int y, std::vector<std::pair<double, double>>& result_V,
    std::vector<std::pair<int, int>>& topo_V);

private:
  std::vector<int> vertices_location_x_flat_;
  std::vector<int> vertices_location_y_flat_;

  std::unordered_map<int, std::vector<std::pair<double, double>>> V_storage_;
  std::unordered_map<int, std::vector<std::pair<int, int>>> topoV_storage_;

  inline int vertexIndex(int x, int y) const { return y * xsize_ + x; }

  bool isInTri(int x1, int y1, int x2, int y2, int x3, int y3, double x, double y);

  inline bool isAnObstacleEdge(std::pair<int, int> prev_pos, std::pair<int, int> next_pos)
  {
    auto topo_prev = locateVertex(prev_pos.first, prev_pos.second);
    auto topo_next = locateVertex(next_pos.first, next_pos.second);
    if (topo_prev.first < 0 || topo_next.first < 0) return false;
    return areConsecutive(topo_prev, topo_next);
  }

  void constructCGALRelated();
  void registerVertices();

  constrained_delaunay_triangulation::CDT cdt_;
  std::unordered_map<long long, int> cdt_table_;
  int cdt_ver_num_;
  std::vector<std::vector<std::pair<int, int>>> facets_;

  inline int locateAdjacentFacet(std::pair<int, int> prev, std::pair<int, int> next) const
  {
    return cdt_table_.at(static_cast<long long>(prev.first + prev.second * xsize_) +
      static_cast<long long>(next.first + next.second * xsize_) * static_cast<long long>(xsize_ * ysize_));
  }
};

}  // namespace raystar
