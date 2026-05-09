#include <raystar/polymap.h>
#include <unordered_set>
#include <stack>
#include <cstring>
#include <cmath>

namespace raystar
{

static std::pair<double, double> findIntersection(
  std::pair<int, int> s, std::pair<int, int> g,
  std::pair<int, int> p, std::pair<int, int> limit)
{
  int a1 = g.second - s.second, b1 = -(g.first - s.first);
  int c1 = -(g.second - s.second) * s.first + (g.first - s.first) * s.second;
  int a2 = limit.second - p.second, b2 = -(limit.first - p.first);
  int c2 = -(limit.second - p.second) * p.first + (limit.first - p.first) * p.second;

  return {
    -1.0 * (c1 * b2 - c2 * b1) / (a1 * b2 - a2 * b1),
    -1.0 * (c1 * a2 - c2 * a1) / (b1 * a2 - b2 * a1)
  };
}

static std::pair<bool, bool> pnpoly(int nvert, double* vertx, double* verty,
  double testx, double testy)
{
  bool rightray = false, leftray = false;
  for (int i = 0, j = nvert - 1; i < nvert; j = i++) {
    double temp = (vertx[j] - vertx[i]) * (testy - verty[i]) / (verty[j] - verty[i]) + vertx[i];
    if (((verty[i] > testy) != (verty[j] > testy)) && (testx < temp))
      rightray = !rightray;
    if (((verty[i] > testy) != (verty[j] > testy)) && (testx > temp))
      leftray = !leftray;
  }
  if (rightray == leftray)
    return {rightray, false};
  return {false, true};
}

bool Polymap::calculateVisibilityRegion(int round_x, int round_y,
  std::vector<std::pair<double, double>>& result_V,
  std::vector<std::pair<int, int>>& topo_V)
{
  std::list<constrained_delaunay_triangulation::BungiuEdge> bd;
  bool open_visibility_region = false;

  int vl_idx = vertices_location_x_flat_[vertexIndex(round_x, round_y)];
  if (vl_idx != -1)
  {
    open_visibility_region = true;
    auto curr = locateVertex(round_x, round_y);
    auto prev = getPrevObs(curr);
    auto next = getNextObs(curr);

    auto the_vertex = prev;
    while (!(the_vertex.first == next.first && the_vertex.second == next.second))
    {
      int loc = locateAdjacentFacet({round_x, round_y}, the_vertex);
      int i;
      for (i = 0; i < 3; ++i) {
        if (!(facets_[loc][i].first == round_x && facets_[loc][i].second == round_y) &&
          !(facets_[loc][i].first == the_vertex.first && facets_[loc][i].second == the_vertex.second))
          break;
      }
      bd.emplace_back(constrained_delaunay_triangulation::BungiuEdge());
      bd.back().prev_pos = the_vertex;
      bd.back().next_pos = facets_[loc][i];
      bd.back().topo_prev = locateVertex(bd.back().prev_pos.first, bd.back().prev_pos.second);
      bd.back().topo_next = locateVertex(bd.back().next_pos.first, bd.back().next_pos.second);
      bd.back().limit_prev = atan2(bd.back().prev_pos.second - round_y, bd.back().prev_pos.first - round_x);
      bd.back().limit_next = atan2(bd.back().next_pos.second - round_y, bd.back().next_pos.first - round_x);
      bd.back().limit_next = bd.back().limit_prev + normalize_angle_positive(bd.back().limit_next - bd.back().limit_prev);
      bd.back().limit_prev_pos = bd.back().prev_pos;
      bd.back().limit_next_pos = bd.back().next_pos;
      bd.back().is_bd = isAnObstacleEdge(bd.back().next_pos, bd.back().prev_pos);
      the_vertex = facets_[loc][i];
    }
  }
  else
  {
    int loc = -1;
    for (size_t fi = 0; fi < facets_.size(); ++fi) {
      if (isInTri(facets_[fi][0].first, facets_[fi][0].second,
        facets_[fi][1].first, facets_[fi][1].second,
        facets_[fi][2].first, facets_[fi][2].second,
        round_x, round_y)) {
        loc = static_cast<int>(fi);
        break;
      }
    }
    for (int i = 0; i < 3; ++i) {
      bd.emplace_back(constrained_delaunay_triangulation::BungiuEdge());
      bd.back().prev_pos = facets_[loc][i];
      bd.back().next_pos = facets_[loc][(i + 1) % 3];
      bd.back().topo_prev = locateVertex(bd.back().prev_pos.first, bd.back().prev_pos.second);
      bd.back().topo_next = locateVertex(bd.back().next_pos.first, bd.back().next_pos.second);
      bd.back().limit_prev = atan2(bd.back().prev_pos.second - round_y, bd.back().prev_pos.first - round_x);
      bd.back().limit_next = atan2(bd.back().next_pos.second - round_y, bd.back().next_pos.first - round_x);
      bd.back().limit_next = bd.back().limit_prev + normalize_angle_positive(bd.back().limit_next - bd.back().limit_prev);
      bd.back().limit_prev_pos = bd.back().prev_pos;
      bd.back().limit_next_pos = bd.back().next_pos;
      bd.back().is_bd = isAnObstacleEdge(bd.back().next_pos, bd.back().prev_pos);
    }
  }

  auto iter = bd.begin();
  double theta;
  int loc;
  while (1)
  {
    if (iter == bd.end()) break;

    if (iter->is_bd) { ++iter; continue; }

    loc = locateAdjacentFacet(iter->next_pos, iter->prev_pos);
    int i;
    for (i = 0; i < 3; ++i) {
      if (!(facets_[loc][i].first == iter->prev_pos.first && facets_[loc][i].second == iter->prev_pos.second) &&
        !(facets_[loc][i].first == iter->next_pos.first && facets_[loc][i].second == iter->next_pos.second))
        break;
    }

    theta = iter->limit_prev + normalize_angle(
      atan2(facets_[loc][i].second - round_y, facets_[loc][i].first - round_x) - iter->limit_prev);

    if (theta < iter->limit_prev)
    {
      if (isAnObstacleEdge(iter->next_pos, facets_[loc][i]))
      {
        auto endpoint_prev = findIntersection(facets_[loc][i], iter->next_pos,
          std::pair<int,int>(round_x, round_y), iter->limit_prev_pos);
        auto endpoint_next = findIntersection(facets_[loc][i], iter->next_pos,
          std::pair<int,int>(round_x, round_y), iter->limit_next_pos);

        auto new_iter = bd.insert(iter, constrained_delaunay_triangulation::BungiuEdge());
        new_iter->prev_pos = iter->limit_prev_pos;
        new_iter->next_pos = endpoint_prev;
        new_iter->topo_prev = locateVertex(new_iter->prev_pos.first, new_iter->prev_pos.second);
        new_iter->topo_next = locateVertex(iter->next_pos.first, iter->next_pos.second);
        new_iter->is_bd = true;

        new_iter = bd.insert(iter, constrained_delaunay_triangulation::BungiuEdge());
        new_iter->prev_pos = endpoint_prev;
        new_iter->next_pos = endpoint_next;
        new_iter->topo_prev = locateVertex(iter->next_pos.first, iter->next_pos.second);
        new_iter->topo_next = locateVertex(iter->next_pos.first, iter->next_pos.second);
        new_iter->is_bd = true;

        if (!(round(iter->next_pos.first) == round(endpoint_next.first) &&
              round(iter->next_pos.second) == round(endpoint_next.second)))
        {
          iter->prev_pos = endpoint_next;
          iter->next_pos = iter->limit_next_pos;
          iter->topo_prev = locateVertex(iter->next_pos.first, iter->next_pos.second);
          iter->topo_next = locateVertex(iter->limit_next_pos.first, iter->limit_next_pos.second);
          iter->is_bd = true;
          iter++;
        } else {
          bd.erase(iter++);
        }
      }
      else
      {
        iter->prev_pos = facets_[loc][i];
        iter->topo_prev = locateVertex(iter->prev_pos.first, iter->prev_pos.second);
      }
    }
    else if (theta == iter->limit_prev)
    {
      auto new_iter = bd.insert(iter, constrained_delaunay_triangulation::BungiuEdge());
      new_iter->prev_pos = iter->limit_prev_pos;
      new_iter->next_pos = facets_[loc][i];
      new_iter->topo_prev = locateVertex(new_iter->prev_pos.first, new_iter->prev_pos.second);
      new_iter->topo_next = locateVertex(new_iter->next_pos.first, new_iter->next_pos.second);
      new_iter->is_bd = true;

      if (isAnObstacleEdge(iter->next_pos, facets_[loc][i]))
      {
        auto endpoint_next = findIntersection(facets_[loc][i], iter->next_pos,
          std::pair<int,int>(round_x, round_y), iter->limit_next_pos);
        new_iter = bd.insert(iter, constrained_delaunay_triangulation::BungiuEdge());
        new_iter->prev_pos = facets_[loc][i];
        new_iter->next_pos = endpoint_next;
        new_iter->topo_prev = locateVertex(new_iter->prev_pos.first, new_iter->prev_pos.second);
        new_iter->topo_next = locateVertex(iter->next_pos.first, iter->next_pos.second);
        new_iter->is_bd = true;

        if (!(round(iter->next_pos.first) == round(endpoint_next.first) &&
              round(iter->next_pos.second) == round(endpoint_next.second)))
        {
          iter->topo_prev = locateVertex(iter->next_pos.first, iter->next_pos.second);
          iter->prev_pos = endpoint_next;
          iter->next_pos = iter->limit_next_pos;
          iter->topo_next = locateVertex(iter->next_pos.first, iter->next_pos.second);
          iter->is_bd = true;
          iter++;
        } else {
          bd.erase(iter++);
        }
      }
      else
      {
        iter->prev_pos = facets_[loc][i];
        iter->topo_prev = locateVertex(iter->prev_pos.first, iter->prev_pos.second);
        iter->limit_prev_pos = facets_[loc][i];
      }
    }
    else if (theta > iter->limit_prev && theta < iter->limit_next)
    {
      bool to_minus = false;
      if (isAnObstacleEdge(facets_[loc][i], iter->prev_pos))
      {
        auto endpoint_prev = findIntersection(iter->prev_pos, facets_[loc][i],
          std::pair<int,int>(round_x, round_y), iter->limit_prev_pos);

        if (!(round(endpoint_prev.first) == round(iter->limit_prev_pos.first) &&
              round(endpoint_prev.second) == round(iter->limit_prev_pos.second)))
        {
          auto ni = bd.insert(iter, constrained_delaunay_triangulation::BungiuEdge());
          ni->prev_pos = iter->limit_prev_pos;
          ni->next_pos = endpoint_prev;
          ni->topo_prev = locateVertex(ni->prev_pos.first, ni->prev_pos.second);
          ni->topo_next = locateVertex(facets_[loc][i].first, facets_[loc][i].second);
          ni->is_bd = true;
        }
        auto ni = bd.insert(iter, constrained_delaunay_triangulation::BungiuEdge());
        ni->prev_pos = endpoint_prev;
        ni->next_pos = facets_[loc][i];
        ni->topo_prev = locateVertex(facets_[loc][i].first, facets_[loc][i].second);
        ni->topo_next = locateVertex(facets_[loc][i].first, facets_[loc][i].second);
        ni->is_bd = true;
      }
      else
      {
        auto ni = bd.insert(iter, constrained_delaunay_triangulation::BungiuEdge());
        ni->prev_pos = iter->prev_pos;
        ni->next_pos = facets_[loc][i];
        ni->topo_prev = locateVertex(facets_[loc][i].first, facets_[loc][i].second);
        ni->topo_next = locateVertex(facets_[loc][i].first, facets_[loc][i].second);
        ni->limit_prev = iter->limit_prev;
        ni->limit_next = theta;
        ni->limit_prev_pos = iter->limit_prev_pos;
        ni->limit_next_pos = facets_[loc][i];
        ni->is_bd = false;
        to_minus = true;
      }

      if (isAnObstacleEdge(iter->next_pos, facets_[loc][i]))
      {
        auto endpoint_next = findIntersection(facets_[loc][i], iter->next_pos,
          std::pair<int,int>(round_x, round_y), iter->limit_next_pos);
        auto ni = bd.insert(iter, constrained_delaunay_triangulation::BungiuEdge());
        ni->prev_pos = facets_[loc][i];
        ni->next_pos = endpoint_next;
        ni->topo_prev = locateVertex(ni->prev_pos.first, ni->prev_pos.second);
        ni->topo_next = locateVertex(iter->next_pos.first, iter->next_pos.second);
        ni->is_bd = true;

        if (!(round(iter->limit_next_pos.first) == round(endpoint_next.first) &&
              round(iter->limit_next_pos.second) == round(endpoint_next.second)))
        {
          iter->prev_pos = endpoint_next;
          iter->topo_prev = locateVertex(iter->next_pos.first, iter->next_pos.second);
          iter->next_pos = iter->limit_next_pos;
          iter->topo_next = locateVertex(iter->next_pos.first, iter->next_pos.second);
          iter->is_bd = true;
          if (to_minus) { iter--; } else { iter++; }
        }
        else
        {
          bd.erase(iter++);
          if (to_minus) { iter--; }
        }
      }
      else
      {
        iter->prev_pos = facets_[loc][i];
        iter->topo_prev = locateVertex(iter->prev_pos.first, iter->prev_pos.second);
        iter->limit_prev = theta;
        iter->limit_prev_pos = facets_[loc][i];
        if (to_minus) iter--;
      }
    }
    else if (theta == iter->limit_next)
    {
      if (isAnObstacleEdge(facets_[loc][i], iter->prev_pos))
      {
        auto endpoint_prev = findIntersection(facets_[loc][i], iter->prev_pos,
          std::pair<int,int>(round_x, round_y), iter->limit_prev_pos);

        if (!(round(endpoint_prev.first) == round(iter->limit_prev_pos.first) &&
              round(endpoint_prev.second) == round(iter->limit_prev_pos.second)))
        {
          auto ni = bd.insert(iter, constrained_delaunay_triangulation::BungiuEdge());
          ni->prev_pos = iter->limit_prev_pos;
          ni->next_pos = endpoint_prev;
          ni->topo_prev = locateVertex(ni->prev_pos.first, ni->prev_pos.second);
          ni->topo_next = locateVertex(facets_[loc][i].first, facets_[loc][i].second);
          ni->is_bd = true;
        }

        auto ni = bd.insert(iter, constrained_delaunay_triangulation::BungiuEdge());
        ni->prev_pos = endpoint_prev;
        ni->next_pos = facets_[loc][i];
        ni->topo_prev = locateVertex(facets_[loc][i].first, facets_[loc][i].second);
        ni->topo_next = locateVertex(facets_[loc][i].first, facets_[loc][i].second);
        ni->is_bd = true;

        iter->prev_pos = facets_[loc][i];
        iter->topo_prev = locateVertex(facets_[loc][i].first, facets_[loc][i].second);
        iter->is_bd = true;
        iter++;
      }
      else
      {
        auto ni = bd.insert(iter, constrained_delaunay_triangulation::BungiuEdge());
        ni->prev_pos = iter->prev_pos;
        ni->next_pos = facets_[loc][i];
        ni->topo_prev = locateVertex(ni->prev_pos.first, ni->prev_pos.second);
        ni->topo_next = locateVertex(ni->next_pos.first, ni->next_pos.second);
        ni->limit_prev = iter->limit_prev;
        ni->limit_next = theta;
        ni->limit_prev_pos = iter->limit_prev_pos;
        ni->limit_next_pos = facets_[loc][i];
        ni->is_bd = false;

        iter->prev_pos = facets_[loc][i];
        iter->topo_prev = locateVertex(facets_[loc][i].first, facets_[loc][i].second);
        iter->is_bd = true;
        iter++;
      }
    }
    else if (theta > iter->limit_next)
    {
      if (isAnObstacleEdge(facets_[loc][i], iter->prev_pos))
      {
        auto endpoint_prev = findIntersection(facets_[loc][i], iter->prev_pos,
          std::pair<int,int>(round_x, round_y), iter->limit_prev_pos);
        auto endpoint_next = findIntersection(facets_[loc][i], iter->prev_pos,
          std::pair<int,int>(round_x, round_y), iter->limit_next_pos);

        if (!(round(endpoint_prev.first) == round(iter->limit_prev_pos.first) &&
              round(endpoint_prev.second) == round(iter->limit_prev_pos.second)))
        {
          auto ni = bd.insert(iter, constrained_delaunay_triangulation::BungiuEdge());
          ni->prev_pos = iter->limit_prev_pos;
          ni->next_pos = endpoint_prev;
          ni->topo_prev = locateVertex(iter->limit_prev_pos.first, iter->limit_prev_pos.second);
          ni->topo_next = locateVertex(facets_[loc][i].first, facets_[loc][i].second);
          ni->is_bd = true;
        }

        auto ni = bd.insert(iter, constrained_delaunay_triangulation::BungiuEdge());
        ni->prev_pos = endpoint_prev;
        ni->next_pos = endpoint_next;
        ni->topo_prev = locateVertex(facets_[loc][i].first, facets_[loc][i].second);
        ni->topo_next = locateVertex(facets_[loc][i].first, facets_[loc][i].second);
        ni->is_bd = true;

        iter->prev_pos = endpoint_next;
        iter->next_pos = iter->limit_next_pos;
        iter->topo_prev = locateVertex(facets_[loc][i].first, facets_[loc][i].second);
        iter->topo_next = locateVertex(iter->limit_next_pos.first, iter->limit_next_pos.second);
        iter->is_bd = true;
        iter++;
      }
      else
      {
        iter->next_pos = facets_[loc][i];
        iter->topo_next = locateVertex(facets_[loc][i].first, facets_[loc][i].second);
      }
    }
    else
    {
      std::cerr << "[raystar] getVisibilityRegion: unexpected branch" << std::endl;
    }
  }

  for (auto it = bd.begin(); it != bd.end(); ++it)
  {
    result_V.emplace_back(it->prev_pos.first, it->prev_pos.second);
    topo_V.emplace_back(it->topo_prev.first, it->topo_prev.second);
  }

  if (open_visibility_region)
  {
    result_V.emplace_back(bd.back().next_pos.first, bd.back().next_pos.second);
    topo_V.emplace_back(bd.back().topo_next.first, bd.back().topo_next.second);
  }

  return !bd.empty();
}

bool Polymap::isInTri(int x1, int y1, int x2, int y2, int x3, int y3, double x, double y)
{
  double vertx[3] = {static_cast<double>(x1), static_cast<double>(x2), static_cast<double>(x3)};
  double verty[3] = {static_cast<double>(y1), static_cast<double>(y2), static_cast<double>(y3)};
  auto b = pnpoly(3, vertx, verty, x, y);
  return b.first || b.second;
}

void Polymap::constructCGALRelated()
{
  for (auto& ob : obs_)
  {
    for (auto it = ob.detailed_ordered_vertices_.begin(); it != ob.detailed_ordered_vertices_.end(); ++it)
    {
      auto next = std::next(it);
      if (next == ob.detailed_ordered_vertices_.end())
        next = ob.detailed_ordered_vertices_.begin();
      cdt_.insert_constraint(
        constrained_delaunay_triangulation::Point(it->first, it->second),
        constrained_delaunay_triangulation::Point(next->first, next->second));
    }
  }

  cdt_.is_valid();

  facets_.clear();
  cdt_ver_num_ = static_cast<int>(cdt_.number_of_vertices());

  int count = 0;
  for (auto fit = cdt_.finite_faces_begin(); fit != cdt_.finite_faces_end(); ++fit)
  {
    facets_.emplace_back(std::vector<std::pair<int,int>>{
      {static_cast<int>(fit->vertex(0)->point().x()), static_cast<int>(fit->vertex(0)->point().y())},
      {static_cast<int>(fit->vertex(1)->point().x()), static_cast<int>(fit->vertex(1)->point().y())},
      {static_cast<int>(fit->vertex(2)->point().x()), static_cast<int>(fit->vertex(2)->point().y())}
    });
    auto& f = facets_.back();
    cdt_table_[static_cast<long long>(f[0].first + f[0].second * xsize_) +
      static_cast<long long>(f[1].first + f[1].second * xsize_) * static_cast<long long>(xsize_ * ysize_)] = count;
    cdt_table_[static_cast<long long>(f[1].first + f[1].second * xsize_) +
      static_cast<long long>(f[2].first + f[2].second * xsize_) * static_cast<long long>(xsize_ * ysize_)] = count;
    cdt_table_[static_cast<long long>(f[2].first + f[2].second * xsize_) +
      static_cast<long long>(f[0].first + f[0].second * xsize_) * static_cast<long long>(xsize_ * ysize_)] = count;
    count++;
  }
}

bool Polymap::getPolyObstacles(int start_x, int start_y, int goal_x, int goal_y)
{
  obs_.clear();
  unsigned int nx = xsize_;
  unsigned int ny = ysize_;

  std::vector<int> mask(nx * ny, 0);

  std::unordered_map<int, int> edges;
  std::stack<int> Q;
  Q.push(start_x + start_y * nx);
  while (!Q.empty())
  {
    int cur = Q.top(); Q.pop();
    int x = cur % nx;
    int y = (cur - x) / nx;
    if (data_[cur] != 0 || mask[x + y * nx] != 0) continue;

    if (cur - 1 >= 0 && data_[cur - 1] != 0)
      edges[cur + cur + nx] = cur;
    if (cur + 1 < static_cast<int>(nx * ny) && data_[cur + 1] != 0)
      edges[cur + 1 + cur + nx + 1] = cur + nx + 1;
    if (cur - static_cast<int>(nx) >= 0 && data_[cur - nx] != 0)
      edges[cur + cur + 1] = cur + 1;
    if (cur + static_cast<int>(nx) < static_cast<int>(nx * ny) && data_[cur + nx] != 0)
      edges[cur + nx + cur + nx + 1] = cur + nx;

    mask[x + y * nx] = 1;
    if (x > 0 && data_[cur - 1] == 0) Q.push(cur - 1);
    if (x < static_cast<int>(nx) - 1 && data_[cur + 1] == 0) Q.push(cur + 1);
    if (y > 0 && data_[cur - nx] == 0) Q.push(cur - nx);
    if (y < static_cast<int>(ny) - 1 && data_[cur + nx] == 0) Q.push(cur + nx);
  }

  if (mask[goal_x + goal_y * nx] == 0)
  {
    std::cout << "[raystar] Polymap construction interrupted: no path to goal" << std::endl;
    return false;
  }

  std::list<std::pair<int, int>> boundary_points;
  while (!edges.empty())
  {
    obs_.emplace_back(Obs());
    boundary_points.clear();

    auto first_iter = edges.begin();
    int key = first_iter->first;
    int value = first_iter->second;

    int prev = value;
    int next = key - value;
    int prev_x = prev % nx;
    int prev_y = (prev - prev_x) / nx;
    int next_x = next % nx;
    int next_y = (next - next_x) / nx;

    boundary_points.emplace_back(prev_x, prev_y);
    boundary_points.emplace_back(next_x, next_y);

    auto biter = boundary_points.end();
    int x, y;
    while (1)
    {
      x = boundary_points.back().first;
      y = boundary_points.back().second;
      int cur = x + y * nx;

      int lb_free = 0, lt_free = 0, rb_free = 0, rt_free = 0;
      if (x > 0 && y > 0 && data_[cur - nx - 1] == 0) lb_free = 1;
      if (x > 0 && data_[cur - 1] == 0) lt_free = 1;
      if (y > 0 && data_[cur - nx] == 0) rb_free = 1;
      if (data_[cur] == 0) rt_free = 1;

      int num = lb_free * 8 + lt_free * 4 + rb_free * 2 + rt_free;

      switch (num)
      {
      case 1: boundary_points.emplace_back(x, y+1); break;
      case 2: boundary_points.emplace_back(x+1, y); break;
      case 3: boundary_points.emplace_back(x, y+1); break;
      case 4: boundary_points.emplace_back(x-1, y); break;
      case 5: boundary_points.emplace_back(x-1, y); break;
      case 6:
        biter = std::prev(boundary_points.end(), 2);
        if ((*biter).first == x && (*biter).second == y - 1)
          boundary_points.emplace_back(x + 1, y);
        else if ((*biter).first == x && (*biter).second == y + 1)
          boundary_points.emplace_back(x - 1, y);
        break;
      case 7: boundary_points.emplace_back(x-1, y); break;
      case 8: boundary_points.emplace_back(x, y-1); break;
      case 9:
        biter = std::prev(boundary_points.end(), 2);
        if ((*biter).first == x + 1 && (*biter).second == y)
          boundary_points.emplace_back(x, y + 1);
        else if ((*biter).first == x - 1 && (*biter).second == y)
          boundary_points.emplace_back(x, y - 1);
        break;
      case 10: boundary_points.emplace_back(x+1, y); break;
      case 11: boundary_points.emplace_back(x, y+1); break;
      case 12: boundary_points.emplace_back(x, y-1); break;
      case 13: boundary_points.emplace_back(x, y-1); break;
      case 14: boundary_points.emplace_back(x+1, y); break;
      }

      if (boundary_points.back().first == boundary_points.front().first &&
        boundary_points.back().second == boundary_points.front().second)
      {
        boundary_points.pop_back();
        break;
      }
    }

    for (auto& bp : boundary_points)
      obs_.back().detailed_ordered_vertices_.emplace_back(bp);

    for (auto it = obs_.back().detailed_ordered_vertices_.begin();
      it != obs_.back().detailed_ordered_vertices_.end(); ++it)
    {
      int curr = it->first + it->second * nx;
      int nxt;
      if (it == obs_.back().detailed_ordered_vertices_.end() - 1)
        nxt = obs_.back().detailed_ordered_vertices_.front().first +
              obs_.back().detailed_ordered_vertices_.front().second * nx;
      else
        nxt = std::next(it)->first + std::next(it)->second * nx;
      edges.erase(curr + nxt);
    }
  }
  return true;
}

void Polymap::getVisibilityRegion(int start_x, int start_y,
  std::vector<std::pair<double, double>>& visibility_region,
  std::vector<std::pair<int, int>>& topo_V)
{
  visibility_region.clear();
  int index = start_x + start_y * xsize_;
  auto it = V_storage_.find(index);
  if (it == V_storage_.end())
  {
    calculateVisibilityRegion(start_x, start_y, visibility_region, topo_V);
    V_storage_[index] = visibility_region;
    topoV_storage_[index] = topo_V;
  }
  else
  {
    visibility_region.assign(it->second.begin(), it->second.end());
    topo_V.assign(topoV_storage_[index].begin(), topoV_storage_[index].end());
  }
}

std::pair<int, int> Polymap::locateVertex(int x, int y) const
{
  int idx = vertexIndex(x, y);
  return {vertices_location_x_flat_[idx], vertices_location_y_flat_[idx]};
}

Polymap::Polymap(const GridMap& grid_map, int start_x, int start_y, int goal_x, int goal_y)
  : xsize_(grid_map.width), ysize_(grid_map.height), data_(grid_map.data)
{
  vertices_location_x_flat_.resize(xsize_ * ysize_, -1);
  vertices_location_y_flat_.resize(xsize_ * ysize_, -1);

  solution_exist_ = getPolyObstacles(start_x, start_y, goal_x, goal_y);
  if (!solution_exist_) return;

  simplifyPolyObstacles(start_x, start_y, goal_x, goal_y);
  registerVertices();
  constructCGALRelated();
}

void Polymap::registerVertices()
{
  for (size_t i = 0; i < obs_.size(); ++i)
  {
    for (size_t j = 0; j < obs_[i].detailed_ordered_vertices_.size(); ++j)
    {
      int x = obs_[i].detailed_ordered_vertices_[j].first;
      int y = obs_[i].detailed_ordered_vertices_[j].second;
      int idx = vertexIndex(x, y);
      vertices_location_x_flat_[idx] = static_cast<int>(i);
      vertices_location_y_flat_[idx] = static_cast<int>(j);
    }
  }
}

void Polymap::simplifyPolyObstacles(int start_x, int start_y, int goal_x, int goal_y)
{
  for (auto iter = obs_.begin(); iter != obs_.end(); ++iter)
  {
    int prev, curr, next;
    bool stable = false;
    curr = 0;
    double prev_dir, next_dir;
    bool simplifable;
    int x1, y1, x2, y2, x3, y3;

    while (1)
    {
      prev = (curr - 1 + static_cast<int>(iter->detailed_ordered_vertices_.size())) %
        static_cast<int>(iter->detailed_ordered_vertices_.size());
      next = (curr + 1) % iter->detailed_ordered_vertices_.size();

      prev_dir = atan2(iter->detailed_ordered_vertices_[curr].second - iter->detailed_ordered_vertices_[prev].second,
        iter->detailed_ordered_vertices_[curr].first - iter->detailed_ordered_vertices_[prev].first);
      next_dir = atan2(iter->detailed_ordered_vertices_[next].second - iter->detailed_ordered_vertices_[curr].second,
        iter->detailed_ordered_vertices_[next].first - iter->detailed_ordered_vertices_[curr].first);

      x1 = iter->detailed_ordered_vertices_[prev].first;
      y1 = iter->detailed_ordered_vertices_[prev].second;
      x2 = iter->detailed_ordered_vertices_[curr].first;
      y2 = iter->detailed_ordered_vertices_[curr].second;
      x3 = iter->detailed_ordered_vertices_[next].first;
      y3 = iter->detailed_ordered_vertices_[next].second;

      if ((x3 - x2) * (y2 - y1) == (x2 - x1) * (y3 - y2))
      {
        simplifable = true;
      }
      else
      {
        double diff_dir = normalize_angle(next_dir - prev_dir);
        simplifable = true;
        if (diff_dir <= 0 || diff_dir > 0.999 * M_PI)
        {
          if (isInTri(x1, y1, x2, y2, x3, y3, start_x, start_y) ||
            isInTri(x1, y1, x2, y2, x3, y3, goal_x, goal_y))
            simplifable = false;

          double testx, testy;
          if (simplifable)
          {
            for (auto iter2 = obs_.begin(); iter2 != obs_.end(); ++iter2)
            {
              for (auto iter3 = iter2->detailed_ordered_vertices_.begin();
                iter3 != iter2->detailed_ordered_vertices_.end(); ++iter3)
              {
                testx = iter3->first;
                testy = iter3->second;
                if (isInTri(x1, y1, x2, y2, x3, y3, testx, testy))
                {
                  if (iter2 - obs_.begin() != iter - obs_.begin())
                    simplifable = false;
                  else
                  {
                    if (iter3 - iter2->detailed_ordered_vertices_.begin() != prev &&
                      iter3 - iter2->detailed_ordered_vertices_.begin() != curr &&
                      iter3 - iter2->detailed_ordered_vertices_.begin() != next)
                      simplifable = false;
                  }
                }
              }
            }
          }
        }
        else
        {
          simplifable = false;
        }
      }

      if (simplifable)
      {
        iter->detailed_ordered_vertices_.erase(iter->detailed_ordered_vertices_.begin() + curr);
        if (curr >= static_cast<int>(iter->detailed_ordered_vertices_.size()))
          curr = static_cast<int>(iter->detailed_ordered_vertices_.size()) - 1;
        stable = false;
        continue;
      }
      else
      {
        if (curr == 0) {
          if (!stable) stable = true;
          else break;
        }
      }

      curr++;
      if (curr >= static_cast<int>(iter->detailed_ordered_vertices_.size()))
        curr = 0;
    }
  }
}

}  // namespace raystar
