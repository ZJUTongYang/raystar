#include <iostream>
#include <vector>
#include <raystar/polymap.h>
#include <algorithm>
#include <unordered_map>
#include <unordered_set>
#include <stack>
#include <raystar/raystar.h>

#include <ros/ros.h>
#include <angles/angles.h>

#include <CGAL/Triangular_expansion_visibility_2.h>

namespace raystar
{
	std::pair<double, double> findIntersection(std::pair<int, int> s, std::pair<int, int> g, std::pair<int, int> p, std::pair<int, int> limit)
	{
		int a1 = g.second - s.second, b1 = -(g.first - s.first), c1 = -(g.second - s.second) * s.first + (g.first - s.first) * s.second;
		int a2 = limit.second - p.second, b2 = -(limit.first - p.first), c2 = -(limit.second - p.second) * p.first + (limit.first - p.first) * p.second;

		std::pair<double, double> result;
		result.first = -1.0 * (c1 * b2 - c2 * b1) / (a1 * b2 - a2 * b1);
		result.second = -1.0 * (c1 * a2 - c2 * a1) / (b1 * a2 - b2 * a1);
		return result;
	}

	std::pair<double, double> findIntersection(std::pair<int, int> s, std::pair<int, int> g, std::pair<double, double> p, double theta)
	{
		double a1 = sin(theta), b1 = -cos(theta), c1 = -sin(theta) * p.first + cos(theta) * p.second;
		double a2 = g.second - s.second, b2 = -(g.first - s.first), c2 = -(g.second - s.second) * s.first + (g.first - s.first) * s.second;

		std::pair<double, double> result;
		result.first = -(c1 * b2 - c2 * b1) / (a1 * b2 - a2 * b1);
		result.second = -(c1 * a2 - c2 * a1) / (b1 * a2 - b2 * a1);
		return result;
	}

bool Polymap::calculateVisibilityRegion(int round_x, int round_y, std::vector<std::pair<double, double> >& result_V, 
	std::vector<std::pair<int, int> >& topo_V)
{
	std::list<constrained_delaunay_triangulation::BungiuEdge> bd;
	bool open_visibility_region = false;
	
	// To create an initial visibility link, we need to check whether the source point is an obstacle vertex
	if (vertices_location_x_[round_x][round_y] != -1)
	{
		// The source point of this visibility region is an obstacle vertex
		open_visibility_region = true;
		std::pair<int, int> curr = locateVertex(round_x, round_y);
		std::pair<int, int> prev = getPrevObs(curr);
		std::pair<int, int> next = getNextObs(curr);

		std::pair<int, int> the_vertex = prev;
		while (!(the_vertex.first == next.first && the_vertex.second == next.second))
		{
			int loc = locateAdjacentFacet(std::pair<double, double>(round_x, round_y), the_vertex);
			int i;
			for (i = 0; i < 3; ++i)
			{
				if (!(facets_[loc][i].first == round_x && facets_[loc][i].second == round_y) &&
					!(facets_[loc][i].first == the_vertex.first && facets_[loc][i].second == the_vertex.second)
					)
					break;
			}
			bd.emplace_back(constrained_delaunay_triangulation::BungiuEdge());
			bd.back().prev_pos = the_vertex;
			bd.back().next_pos = facets_[loc][i];
			bd.back().topo_prev = locateVertex(bd.back().prev_pos.first, bd.back().prev_pos.second);
			bd.back().topo_next = locateVertex(bd.back().next_pos.first, bd.back().next_pos.second);
			bd.back().limit_prev = atan2(bd.back().prev_pos.second - round_y, bd.back().prev_pos.first - round_x);
			bd.back().limit_next = atan2(bd.back().next_pos.second - round_y, bd.back().next_pos.first - round_x);
			bd.back().limit_next = bd.back().limit_prev + angles::normalize_angle_positive(bd.back().limit_next - bd.back().limit_prev);
			bd.back().limit_prev_pos = bd.back().prev_pos;
			bd.back().limit_next_pos = bd.back().next_pos;
			bd.back().is_bd = isAnObstacleEdge(bd.back().next_pos, bd.back().prev_pos);
			the_vertex = facets_[loc][i];
		}
	}
	else
	{
		// The source point is the starting location. It must be within a triangular facet.
		int loc = -1;
		for (auto iter = facets_.begin(); iter != facets_.end(); ++iter)
		{
			if (isInTri((*iter)[0].first, (*iter)[0].second, (*iter)[1].first, (*iter)[1].second, (*iter)[2].first, (*iter)[2].second, round_x, round_y))
			{
				loc = iter - facets_.begin();
				break;
			}
		}
		for (unsigned int i = 0; i < 3; ++i)
		{
			bd.emplace_back(constrained_delaunay_triangulation::BungiuEdge());
			bd.back().prev_pos = facets_[loc][i];
			bd.back().next_pos = facets_[loc][(i + 1) % 3];
			bd.back().topo_prev = locateVertex(bd.back().prev_pos.first, bd.back().prev_pos.second);
			bd.back().topo_next = locateVertex(bd.back().next_pos.first, bd.back().next_pos.second);
			bd.back().limit_prev = atan2(bd.back().prev_pos.second - round_y, bd.back().prev_pos.first - round_x);
			bd.back().limit_next = atan2(bd.back().next_pos.second - round_y, bd.back().next_pos.first - round_x);
			bd.back().limit_next = bd.back().limit_prev + angles::normalize_angle_positive(bd.back().limit_next - bd.back().limit_prev);
			bd.back().limit_prev_pos = bd.back().prev_pos;
			bd.back().limit_next_pos = bd.back().next_pos;
			bd.back().is_bd = isAnObstacleEdge(bd.back().next_pos, bd.back().prev_pos);
		}
	}

	// We iteratively expand the link of visibility boundary.
	auto iter = bd.begin();
	double theta;
	int loc;
	while (1)
	{
		if (iter == bd.end())
			break;

		if (iter->is_bd)
		{
			++iter;
			continue;
		}

		// We find the adjacent facet
		loc = locateAdjacentFacet(iter->next_pos, iter->prev_pos);
		int i;
		for (i = 0; i < 3; ++i)
		{
			if (!(facets_[loc][i].first == iter->prev_pos.first && facets_[loc][i].second == iter->prev_pos.second) &&
				!(facets_[loc][i].first == iter->next_pos.first && facets_[loc][i].second == iter->next_pos.second)
				)
				break;
		}

		// We calculate the new angle
		theta = iter->limit_prev + angles::normalize_angle(atan2(facets_[loc][i].second - round_y, facets_[loc][i].first - round_x) - iter->limit_prev);

		if (theta < iter->limit_prev)
		{
			if (isAnObstacleEdge(iter->next_pos, facets_[loc][i]))
			{
				std::pair<double, double> endpoint_prev = findIntersection(facets_[loc][i], iter->next_pos, std::pair<int, int>(round_x, round_y), iter->limit_prev_pos);
				std::pair<double, double> endpoint_next = findIntersection(facets_[loc][i], iter->next_pos, std::pair<int, int>(round_x, round_y), iter->limit_next_pos);

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

				if (!(round(iter->next_pos.first) == round(endpoint_next.first) && round(iter->next_pos.second) == round(endpoint_next.second)))
				{
					iter->prev_pos = endpoint_next;
					iter->next_pos = iter->limit_next_pos;
					iter->topo_prev = locateVertex(iter->next_pos.first, iter->next_pos.second);
					iter->topo_next = locateVertex(iter->limit_next_pos.first, iter->limit_next_pos.second);
					iter->is_bd = true;
					iter++;
				}
				else
				{
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
				std::pair<double, double> endpoint_next = findIntersection(facets_[loc][i], iter->next_pos, std::pair<int, int>(round_x, round_y), iter->limit_next_pos);
				auto new_iter = bd.insert(iter, constrained_delaunay_triangulation::BungiuEdge());
				new_iter->prev_pos = facets_[loc][i];
				new_iter->next_pos = endpoint_next;
				new_iter->topo_prev = locateVertex(new_iter->prev_pos.first, new_iter->prev_pos.second);
				new_iter->topo_next = locateVertex(iter->next_pos.first, iter->next_pos.second);
				new_iter->is_bd = true;

				if (!(round(iter->next_pos.first) == round(endpoint_next.first) && round(iter->next_pos.second) == round(endpoint_next.second)))
				{
					iter->topo_prev = locateVertex(iter->next_pos.first, iter->next_pos.second);
					iter->prev_pos = endpoint_next;
					iter->next_pos = iter->limit_next_pos;
					iter->topo_next = locateVertex(iter->next_pos.first, iter->next_pos.second);
					iter->is_bd = true;
					iter++;
				}
				else
				{
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
				std::pair<double, double> endpoint_prev = findIntersection(iter->prev_pos, facets_[loc][i], std::pair<int, int>(round_x, round_y), iter->limit_prev_pos);

				if (!(round(endpoint_prev.first) == round(iter->limit_prev_pos.first) && round(endpoint_prev.second) == round(iter->limit_prev_pos.second)))
				{
					auto new_iter = bd.insert(iter, constrained_delaunay_triangulation::BungiuEdge());
					new_iter->prev_pos = iter->limit_prev_pos;
					new_iter->next_pos = endpoint_prev;
					new_iter->topo_prev = locateVertex(new_iter->prev_pos.first, new_iter->prev_pos.second);
					new_iter->topo_next = locateVertex(facets_[loc][i].first, facets_[loc][i].second);
					new_iter->is_bd = true;
				}
				auto new_iter = bd.insert(iter, constrained_delaunay_triangulation::BungiuEdge());
				new_iter->prev_pos = endpoint_prev;
				new_iter->next_pos = facets_[loc][i];
				new_iter->topo_prev = locateVertex(facets_[loc][i].first, facets_[loc][i].second);
				new_iter->topo_next = locateVertex(facets_[loc][i].first, facets_[loc][i].second);
				new_iter->is_bd = true;
			}
			else
			{
				auto new_iter = bd.insert(iter, constrained_delaunay_triangulation::BungiuEdge());
				new_iter->prev_pos = iter->prev_pos;
				new_iter->next_pos = facets_[loc][i];
				new_iter->topo_prev = locateVertex(facets_[loc][i].first, facets_[loc][i].second);
				new_iter->topo_next = locateVertex(facets_[loc][i].first, facets_[loc][i].second);
				new_iter->limit_prev = iter->limit_prev;
				new_iter->limit_next = theta;
				new_iter->limit_prev_pos = iter->limit_prev_pos;
				new_iter->limit_next_pos = facets_[loc][i];
				new_iter->is_bd = false;
				to_minus = true;
			}

			if (isAnObstacleEdge(iter->next_pos, facets_[loc][i]))
			{
				std::pair<double, double> endpoint_next = findIntersection(facets_[loc][i], iter->next_pos, std::pair<int, int>(round_x, round_y), iter->limit_next_pos);
				auto new_iter = bd.insert(iter, constrained_delaunay_triangulation::BungiuEdge());
				new_iter->prev_pos = facets_[loc][i];
				new_iter->next_pos = endpoint_next;
				new_iter->topo_prev = locateVertex(new_iter->prev_pos.first, new_iter->prev_pos.second);
				new_iter->topo_next = locateVertex(iter->next_pos.first, iter->next_pos.second);
				new_iter->is_bd = true;

				if (!(round(iter->limit_next_pos.first) == round(endpoint_next.first) && round(iter->limit_next_pos.second) == round(endpoint_next.second)))
				{
					iter->prev_pos = endpoint_next;
					iter->topo_prev = locateVertex(iter->next_pos.first, iter->next_pos.second);
					iter->next_pos = iter->limit_next_pos;
					iter->topo_next = locateVertex(iter->next_pos.first, iter->next_pos.second);
					iter->is_bd = true;
					if (to_minus)
					{
						iter = --new_iter;
					}
					else
					{
						iter++;
					}
				}
				else
				{
					bd.erase(iter++);
					if (to_minus)
					{
						iter = --new_iter;
					}
				}
			}
			else
			{
				iter->prev_pos = facets_[loc][i];
				iter->topo_prev = locateVertex(iter->prev_pos.first, iter->prev_pos.second);
				iter->limit_prev = theta;
				iter->limit_prev_pos = facets_[loc][i];
				if (to_minus)
					iter--;
			}

		}
		else if (theta == iter->limit_next)
		{
			if (isAnObstacleEdge(facets_[loc][i], iter->prev_pos))
			{
				std::pair<double, double> endpoint_prev = findIntersection(facets_[loc][i], iter->prev_pos, std::pair<int, int>(round_x, round_y), iter->limit_prev_pos);

				if (!(round(endpoint_prev.first) == round(iter->limit_prev_pos.first) && round(endpoint_prev.second) == round(iter->limit_prev_pos.second)))
				{
					auto new_iter = bd.insert(iter, constrained_delaunay_triangulation::BungiuEdge());
					new_iter->prev_pos = iter->limit_prev_pos;
					new_iter->next_pos = endpoint_prev;
					new_iter->topo_prev = locateVertex(new_iter->prev_pos.first, new_iter->prev_pos.second);
					new_iter->topo_next = locateVertex(facets_[loc][i].first, facets_[loc][i].second);
					new_iter->is_bd = true;
				}

				auto new_iter = bd.insert(iter, constrained_delaunay_triangulation::BungiuEdge());
				new_iter->prev_pos = endpoint_prev;
				new_iter->next_pos = facets_[loc][i];
				new_iter->topo_prev = locateVertex(facets_[loc][i].first, facets_[loc][i].second);
				new_iter->topo_next = locateVertex(facets_[loc][i].first, facets_[loc][i].second);
				new_iter->is_bd = true;

				iter->prev_pos = facets_[loc][i];
				iter->topo_prev = locateVertex(facets_[loc][i].first, facets_[loc][i].second);
				iter->is_bd = true;
				iter++;
			}
			else
			{
				auto new_iter = bd.insert(iter, constrained_delaunay_triangulation::BungiuEdge());
				new_iter->prev_pos = iter->prev_pos;
				new_iter->next_pos = facets_[loc][i];
				new_iter->topo_prev = locateVertex(new_iter->prev_pos.first, new_iter->prev_pos.second);
				new_iter->topo_next = locateVertex(new_iter->next_pos.first, new_iter->next_pos.second);
				new_iter->limit_prev = iter->limit_prev;
				new_iter->limit_next = theta;
				new_iter->limit_prev_pos = iter->limit_prev_pos;
				new_iter->limit_next_pos = facets_[loc][i];
				new_iter->is_bd = false;

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
				std::pair<double, double> endpoint_prev = findIntersection(facets_[loc][i], iter->prev_pos, std::pair<int, int>(round_x, round_y), iter->limit_prev_pos);
				std::pair<double, double> endpoint_next = findIntersection(facets_[loc][i], iter->prev_pos, std::pair<int, int>(round_x, round_y), iter->limit_next_pos);

				if (!(round(endpoint_prev.first) == round(iter->limit_prev_pos.first) && round(endpoint_prev.second) == round(iter->limit_prev_pos.second)))
				{
					auto new_iter = bd.insert(iter, constrained_delaunay_triangulation::BungiuEdge());
					new_iter->prev_pos = iter->limit_prev_pos;
					new_iter->next_pos = endpoint_prev;
					new_iter->topo_prev = locateVertex(iter->limit_prev_pos.first, iter->limit_prev_pos.second);
					new_iter->topo_next = locateVertex(facets_[loc][i].first, facets_[loc][i].second);
					new_iter->is_bd = true;
				}

				auto new_iter = bd.insert(iter, constrained_delaunay_triangulation::BungiuEdge());
				new_iter->prev_pos = endpoint_prev;
				new_iter->next_pos = endpoint_next;
				new_iter->topo_prev = locateVertex(facets_[loc][i].first, facets_[loc][i].second);
				new_iter->topo_next = locateVertex(facets_[loc][i].first, facets_[loc][i].second);
				new_iter->is_bd = true;

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
			std::cout << "getVisibilityRegion: We encounter a bug" << std::endl;
		}
	}// while

	for (auto iter = bd.begin(); iter != bd.end(); ++iter)
	{
		result_V.emplace_back(std::pair<double, double>(iter->prev_pos.first, iter->prev_pos.second));
		topo_V.emplace_back(std::pair<int, int>(iter->topo_prev.first, iter->topo_prev.second));
	}

	if (open_visibility_region)
	{
		result_V.emplace_back(std::pair<double, double>(bd.back().next_pos.first, bd.back().next_pos.second));
		topo_V.emplace_back(std::pair<int, int>(bd.back().topo_next.first, bd.back().topo_next.second));
	}

	return !bd.empty();
}

bool Polymap::isInTri(int x1, int y1, int x2, int y2, int x3, int y3, double x, double y)
{
	double vertx[3] = { x1, x2, x3 };
	double verty[3] = { y1, y2, y3 };
	std::pair<bool, bool> b = pnpoly(3, vertx, verty, x, y);
	return b.first || b.second;
}

void Polymap::constructCGALRelated()
{
	for (auto iter = obs_.begin(); iter != obs_.end(); ++iter)
	{
		for (auto iter2 = iter->detailed_ordered_vertices_.begin(); iter2 != iter->detailed_ordered_vertices_.end(); ++iter2)
		{
			auto next = std::next(iter2);
			if (next == iter->detailed_ordered_vertices_.end())
			{
				next = iter->detailed_ordered_vertices_.begin();
			}
			cdt_.insert_constraint(constrained_delaunay_triangulation::Point(iter2->first, iter2->second), constrained_delaunay_triangulation::Point(next->first, next->second));
		}
	}

	cdt_.is_valid();

	std::vector<std::vector<std::pair<int, int> > > facets;
	facets.clear();
	std::vector<std::pair<int, int> > temp;

	cdt_ver_num_ = cdt_.number_of_vertices();

	int x, y;
	int count = 0;
	for (constrained_delaunay_triangulation::CDT::Finite_faces_iterator fit = cdt_.finite_faces_begin(); fit != cdt_.finite_faces_end(); ++fit) 
	{
		facets.emplace_back(temp);
		facets.back().emplace_back(std::pair<int, int>(fit->vertex(0)->point().x(), fit->vertex(0)->point().y()));
		facets.back().emplace_back(std::pair<int, int>(fit->vertex(1)->point().x(), fit->vertex(1)->point().y()));
		facets.back().emplace_back(std::pair<int, int>(fit->vertex(2)->point().x(), fit->vertex(2)->point().y()));
		cdt_table_[(long long)(fit->vertex(0)->point().x() + fit->vertex(0)->point().y() * xsize_) + (long long)(fit->vertex(1)->point().x() + fit->vertex(1)->point().y() * xsize_) * (long long)(xsize_*ysize_)] = count;
		cdt_table_[(long long)(fit->vertex(1)->point().x() + fit->vertex(1)->point().y() * xsize_) + (long long)(fit->vertex(2)->point().x() + fit->vertex(2)->point().y() * xsize_) * (long long)(xsize_ * ysize_)] = count;
		cdt_table_[(long long)(fit->vertex(2)->point().x() + fit->vertex(2)->point().y() * xsize_) + (long long)(fit->vertex(0)->point().x() + fit->vertex(0)->point().y() * xsize_) * (long long)(xsize_ * ysize_)] = count;

		count++;
	}
	facets_ = std::move(facets);
}

bool Polymap::getPolyObstacles(int start_x, int start_y, int goal_x, int goal_y)
{
	obs_.clear();
	unsigned int nx = xsize_;
	unsigned int ny = ysize_;

	int** mask = new int*[nx];
	for (unsigned int i = 0; i < nx; ++i)
	{
		mask[i] = new int[ny];
		memset(mask[i], 0, ny * sizeof(int));
	}

	// (1) We find the collision-free environment. Altogether, we identify obstacle boundaries
	std::unordered_map<int, int> edges; // <source+target, source>
	std::stack<int> Q;	
	Q.emplace(start_x + start_y * nx);
	while(!Q.empty())
	{
		const auto cur = Q.top();
		int x = cur % nx;
		int y = (cur - x) / nx;
		Q.pop();
		if (data_[cur] != 0 || mask[x][y] != 0)
		{
			continue;
		}
		else
		{
			if (data_[cur-1] != 0)
				edges[cur + cur + nx] = cur; // nx
			if (data_[cur + 1] != 0)
				edges[cur + 1 + cur + nx + 1] = cur + nx + 1; // -nx
			if (data_[cur - nx] != 0)
				edges[cur + cur + 1] = cur + 1; // -1
			if (data_[cur +nx] != 0)
				edges[cur + nx + cur + nx + 1] = cur + nx; // 1

			mask[x][y] = 1;
			if (data_[cur - 1] == 0)
				Q.push(cur - 1);
			if (data_[cur + 1] == 0)
				Q.push(cur + 1);
			if (data_[cur - nx] == 0)
				Q.push(cur - nx);
			if (data_[cur + nx] == 0)
				Q.push(cur + nx);
		}
	}

	// (2) If the goal position is not reachable, we don't need to do planning actually
	if (mask[goal_x][goal_y] == 0)
	{
		std::cout << "Polymap construction interrupted as there will be no result path" << std::endl;
		for (unsigned int i = 0; i < nx; ++i)
		{
			delete[] mask[i];
		}
		delete[] mask;
		return false;
	}

	for (unsigned int i = 0; i < nx; ++i)
	{
		delete[] mask[i];
	}
	delete[] mask;

	// (3) We identify the boundary of each obstacle using Bug algorithm
	std::list<std::pair<int, int> > boundary_points;
	while(!edges.empty())
	{
		obs_.emplace_back(Obs());

		boundary_points.clear();

		auto first_iter = edges.begin();
		int key = first_iter->first;
		int value = first_iter->second;

		int dir = key - 2*value; // next - prev
		int prev = value;
		int next = key - value;
		int prev_x = prev % nx;
		int prev_y = (prev - prev_x) / nx; // x, y are calculated based on prev
		int next_x = next % nx;
		int next_y = (next - next_x) / nx;
		
		boundary_points.emplace_back(std::pair<int, int>(prev_x, prev_y));
		boundary_points.emplace_back(std::pair<int, int>(next_x, next_y));

		auto iter = boundary_points.end();
		int x, y;
		while(1)
		{
			x = boundary_points.back().first;
			y = boundary_points.back().second;
			int cur = x + y * nx;

			int lb_free = 0, lt_free = 0, rb_free = 0, rt_free = 0;

			if(x > 0 && y > 0 && data_[cur - nx - 1] == 0)
				lb_free = 1; // left bottom
			if(x > 0 && data_[cur - 1] == 0)
				lt_free = 1; // left top
			if(y > 0 && data_[cur - nx] == 0)
				rb_free = 1; // right bottom
			if(data_[cur] == 0)
				rt_free = 1; // right top

			// a binary number representing [lb, lt, rb, rt] 
			int num = lb_free * 8 + lt_free * 4 + rb_free * 2 + rt_free;

			// position map
			//  2 | 4
			// -------
			//  1 | 3

			switch (num)
			{
			case 1:
				boundary_points.emplace_back(std::pair<int, int>(x, y+1));
				break;
			case 2: 
				boundary_points.emplace_back(std::pair<int, int>(x+1, y));
				break;
			case 3: 
				boundary_points.emplace_back(std::pair<int, int>(x, y+1));
				break;
			case 4: 
				boundary_points.emplace_back(std::pair<int, int>(x-1, y));
				break;
			case 5: 
				boundary_points.emplace_back(std::pair<int, int>(x-1, y));
				break;
			case 6: 
				iter = std::prev(boundary_points.end(), 2);
				if ((*iter).first == x && (*iter).second == y - 1) { 
					boundary_points.emplace_back(std::pair<int, int>(x + 1, y));
				}
				else if ((*iter).first == x && (*iter).second == y + 1) { 
					boundary_points.emplace_back(std::pair<int, int>(x - 1, y));
				}
				break;
			case 7: 
				boundary_points.emplace_back(std::pair<int, int>(x-1, y));
				break;
			case 8: 
				boundary_points.emplace_back(std::pair<int, int>(x, y-1));
				break;
			case 9:
				iter = std::prev(boundary_points.end(), 2);
				if ((*iter).first == x + 1 && (*iter).second == y) { 
					boundary_points.emplace_back(std::pair<int, int>(x, y + 1));
				}
				else if ((*iter).first == x - 1 && (*iter).second == y) { 
					boundary_points.emplace_back(std::pair<int, int>(x, y - 1));
				}
				break;
			case 10: 
				boundary_points.emplace_back(std::pair<int, int>(x+1, y));
				break;
			case 11: 
				boundary_points.emplace_back(std::pair<int, int>(x, y+1));
				break;
			case 12: 
				boundary_points.emplace_back(std::pair<int, int>(x, y-1));
				break;
			case 13: 
				boundary_points.emplace_back(std::pair<int, int>(x, y-1));
				break;
			case 14: 
				boundary_points.emplace_back(std::pair<int, int>(x+1, y));
				break;
			}

			if (boundary_points.back().first == boundary_points.front().first && 
				boundary_points.back().second == boundary_points.front().second)
			{
				boundary_points.pop_back();
				break;
			}
		}

		for(auto iter = boundary_points.begin(); iter != boundary_points.end(); ++iter)
			obs_.back().detailed_ordered_vertices_.emplace_back(*iter);

		// We remove the used edges
		for(auto iter = obs_.back().detailed_ordered_vertices_.begin(); iter != obs_.back().detailed_ordered_vertices_.end(); ++iter)
		{
			int curr = iter->first + iter->second*nx;
			int next;
			if(iter == obs_.back().detailed_ordered_vertices_.end()-1)
			{
				next = obs_.back().detailed_ordered_vertices_.front().first + obs_.back().detailed_ordered_vertices_.front().second * nx;
			}
			else
			{
				next = std::next(iter)->first + std::next(iter)->second * nx;
			}

			edges.erase(curr+next);
		}
	}
	return true;
}

void Polymap::getVisibilityRegion(int start_x, int start_y,
	std::vector<std::pair<double, double> >& visibility_region,
	std::vector<std::pair<int, int> >& topo_V)
{
	visibility_region.clear();
	
	//We try reusing the visibility regions
	int index = start_x + start_y * xsize_;
	auto iter = V_storage_.find(index);
	if (iter == V_storage_.end())
	{
		calculateVisibilityRegion(start_x, start_y, visibility_region, topo_V);
		V_storage_[index] = visibility_region;
		topoV_storage_[index] = topo_V;
	}
	else
	{
		visibility_region.assign(iter->second.begin(), iter->second.end());
		auto& topo_iter = topoV_storage_[index];
		topo_V.assign(topo_iter.begin(), topo_iter.end());
	}
}

std::pair<int, int> Polymap::locateVertex(int x, int y) const
{
	return std::pair<int, int>(vertices_location_x_[x][y], vertices_location_y_[x][y]);
}

Polymap::Polymap(costmap_2d::Costmap2D* costmap, int start_x, int start_y, int goal_x, int goal_y)
{
    data_ = costmap->getCharMap();
	this->xsize_ = costmap->getSizeInCellsX();
	this->ysize_ = costmap->getSizeInCellsY();

	vertices_location_x_ = new int* [xsize_];
	vertices_location_y_ = new int* [xsize_];
	for (unsigned int i = 0; i < xsize_; ++i)
	{
		vertices_location_x_[i] = new int[ysize_];
		vertices_location_y_[i] = new int[ysize_];
		memset(vertices_location_x_[i], -1, sizeof(int) * ysize_);
		memset(vertices_location_y_[i], -1, sizeof(int) * ysize_);
	}

	// When we identify the collision-free space from the robot's current location, we will know the existence of solutions	
	solution_exist_ = getPolyObstacles(start_x, start_y, goal_x, goal_y);

	if(!solution_exist_)
		return ;

	simplifyPolyObstacles(start_x, start_y, goal_x, goal_y);

	registerVertices();

	constructCGALRelated();
}

Polymap::~Polymap()
{
	for (unsigned int i = 0; i < xsize_; ++i)
	{
		delete[] vertices_location_x_[i];
		delete[] vertices_location_y_[i];
	}
	delete[] vertices_location_x_;
	delete[] vertices_location_y_;
}

void Polymap::registerVertices()
{
	for(unsigned int i = 0; i < obs_.size(); ++i)
	{
		for(unsigned int j = 0; j < obs_[i].detailed_ordered_vertices_.size(); ++j)
		{
			int x = obs_[i].detailed_ordered_vertices_[j].first;
			int y = obs_[i].detailed_ordered_vertices_[j].second;
			vertices_location_x_[x][y] = i;
			vertices_location_y_[x][y] = j;
		}
	}
}

void Polymap::simplifyPolyObstacles(int start_x, int start_y, int goal_x, int goal_y)
{
	for(auto iter = obs_.begin(); iter != obs_.end(); ++iter)
	{
		// YT: we simplify obstacle vertices
		int prev, curr, next;
		bool stable = false;
		curr = 0;
		double prev_dir, next_dir, diff_dir;
		bool simplifable;
		int x1, y1, x2, y2, x3, y3;

		while(1)
		{
			prev = (curr - 1 + iter->detailed_ordered_vertices_.size()) % iter->detailed_ordered_vertices_.size();

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
				// If three obstacle vertices are parallel, we directly simplify
				simplifable = true;
			}
			else
			{
				diff_dir = angles::normalize_angle(next_dir - prev_dir);

				simplifable = true;
				if (diff_dir <= 0 || diff_dir > 0.999 * M_PI)
				{
					// Start point and the goal point must not be in the inflated area
					if (isInTri(x1, y1, x2, y2, x3, y3, start_x, start_y) || isInTri(x1, y1, x2, y2, x3, y3, goal_x, goal_y))
						simplifable = false;

					double testx, testy;
					if (simplifable)
					{
						// Check that no other obstacle vertex in the inflated area
						for (auto iter2 = obs_.begin(); iter2 != obs_.end(); ++iter2)
						{
							for (auto iter3 = iter2->detailed_ordered_vertices_.begin(); iter3 != iter2->detailed_ordered_vertices_.end(); ++iter3)
							{
								testx = iter3->first;
								testy = iter3->second;
								if (isInTri(x1, y1, x2, y2, x3, y3, testx, testy))
								{
									// The triangle vertices themselves do not count
									if (iter2 - obs_.begin() != iter - obs_.begin())
										simplifable = false;
									else
									{
										if (iter3 - iter2->detailed_ordered_vertices_.begin() != prev &&
											iter3 - iter2->detailed_ordered_vertices_.begin() != curr &&
											iter3 - iter2->detailed_ordered_vertices_.begin() != next)
										{
											simplifable = false;
										}
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

			if(simplifable)
			{
				iter->detailed_ordered_vertices_.erase(iter->detailed_ordered_vertices_.begin()+curr);
				if(curr >= iter->detailed_ordered_vertices_.size())
					curr = iter->detailed_ordered_vertices_.size()-1;
				stable = false;
				continue;
			}
			else
			{
				if(curr == 0)
				{
					if(!stable)
						stable = true;
					else
						break;
				}
			}

			curr++;
			if(curr >= iter->detailed_ordered_vertices_.size())
			{
				curr = 0;
			}
		}
	}
}
}
