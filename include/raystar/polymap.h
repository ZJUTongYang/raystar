#pragma once

#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <raystar/raystar.h>
#include <costmap_2d/costmap_2d.h>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Triangulation_vertex_base_with_info_2.h>
#include <CGAL/Triangulation_face_base_with_info_2.h>

namespace raystar
{

/* When constructing visibility region, we adopt this paper: 
   F. Bungiu, et al., Efficient computation of visibility polygons, CoRR, abs/1403.3905, 2014
*/ 
namespace constrained_delaunay_triangulation
{
	// Although we adopt CGAL, the implementation of Triangular_expansion_visibility_2 is not usable when the source point is an obstacle vertex. 
	// Therefore, we adopt Constrained_Delaunay_triangulation_2 from CGAL, and implement 2D visibility region by ourselves.
	typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
	typedef CGAL::Triangulation_vertex_base_2<K> Vb;
	typedef CGAL::Triangulation_face_base_with_info_2<int, K> Fb;
	typedef CGAL::Constrained_Delaunay_triangulation_2<K, CGAL::Default, CGAL::No_constraint_intersection_tag> CDT;
	typedef CDT::Point Point;
	typedef CDT::Vertex_handle Vertex_handle;
	typedef CDT::Face_handle Face_handle;
	
	// Refer to Fig. 1 of Bungiu's paper for the meaning of these parameters
	struct BungiuEdge
	{
		// the "a"
		std::pair<double, double> prev_pos;
		// the "b"
		std::pair<double, double> next_pos;
		// the index of the vertex "a"
		std::pair<int, int> topo_prev;
		// the index of the vertex "b"
		std::pair<int, int> topo_next;
		// the orientation of the dashed blue line segment at r
		double limit_prev;
		// the orientation of the dashed blue line segment at l
		double limit_next;
		// the "r"
		std::pair<int, int> limit_prev_pos;
		// the "l"
		std::pair<int, int> limit_next_pos;
		// the flag indicating whether "ab" is an obstacle edge
		bool is_bd;
	};
}

class Obs
{
public:
	// If you want to calculate homotopy invariants such as H-signature, you need to set rep_point_ to be one of the obstacle vertices. 
	// This is not necessary for Ray*.
	//std::pair<int, int> rep_point_;

	// Location of obstacle vertices in order
	std::vector<std::pair<int, int> > detailed_ordered_vertices_;
};

class Polymap
{
public:
	// interface to the costmap_2d in ROS
	int xsize_;
	int ysize_;
    unsigned char* data_;

	std::vector<Obs> obs_;

	// After we identify the obstacles in the environment, whether a result path exists can be determined. 
	// We will cancel useless calculation if no path will be found.  
	bool solution_exist_;
    
	Polymap(costmap_2d::Costmap2D* costmap, 
		int start_x, int start_y, int goal_x, int goal_y);
	~Polymap();

	//std::vector<int> calSwing(const std::vector<int>& oldH, int oldx, int oldy, int newx, int newy) const;

	// Transform the grid-based map into a polygonal map
	bool getPolyObstacles(int start_x, int start_y, int goal_x, int goal_y);

	/* The source point (start_x, start_y) must be an integer point,
	   while the visibility region may contain float-point numbers.
	   topo_V is the index of the visibility region vertices
	*/
	void getVisibilityRegion(int start_x, int start_y, 
		std::vector<std::pair<double, double> >& visibility_region,
		std::vector<std::pair<int, int> >& topo_V);

	// If an obstacle is locally concave, and the concaved area does not contain any obstacle vertices, nor the start / goal point, then we can "convexize" it. 
	void simplifyPolyObstacles(int start_x, int start_y, int goal_x, int goal_y);

	// get the index of an obstacle vertex
	std::pair<int, int> locateVertex(int x, int y) const;

	// check whether two obstacle vertices form an edge of an obstacle
	inline bool areConsecutive(const std::pair<int, int>& prev, const std::pair<int, int>& next) const
	{
		return (prev.first == next.first) && 
			(next.second - prev.second + obs_[prev.first].detailed_ordered_vertices_.size()) % obs_[prev.first].detailed_ordered_vertices_.size() == 1;
	}

	inline std::pair<int, int> getPrevObs(const std::pair<int, int>& curr) const
	{
		return obs_[curr.first].detailed_ordered_vertices_[(curr.second - 1 + obs_[curr.first].detailed_ordered_vertices_.size()) % obs_[curr.first].detailed_ordered_vertices_.size()];
	}
	inline std::pair<int, int> getNextObs(const std::pair<int, int>& curr) const
	{
		return obs_[curr.first].detailed_ordered_vertices_[(curr.second + 1) % obs_[curr.first].detailed_ordered_vertices_.size()];
	}

	bool calculateVisibilityRegion(int x, int y, std::vector<std::pair<double, double> >& result_V,
		std::vector<std::pair<int, int> >& topo_V);

private: 
	// Storage of the index of obstacle vertices, for fast query
	int** vertices_location_x_;
	int** vertices_location_y_;

	// The algorithm may calculate the visibility region on the same source point for multiple times. 
	// Hence we store the V-regions, for fast multi-query. 
	std::unordered_map<int, std::vector<std::pair<double, double> > > V_storage_;
	std::unordered_map<int, std::vector<std::pair<int, int> > > topoV_storage_;

	// Check whether the query point is on the boundary of or within the triangle
	bool isInTri(int x1, int y1, int x2, int y2, int x3, int y3, double x, double y);

	// For two obstacle vertices (not their indices), check whether they form an obstacle edge
	inline bool isAnObstacleEdge(std::pair<int, int> prev_pos, std::pair<int, int> next_pos)
	{
		std::pair<int, int> topo_prev(vertices_location_x_[prev_pos.first][prev_pos.second], vertices_location_y_[prev_pos.first][prev_pos.second]);
		std::pair<int, int> topo_next(vertices_location_x_[next_pos.first][next_pos.second], vertices_location_y_[next_pos.first][next_pos.second]);
		return areConsecutive(topo_prev, topo_next);
	}

	void constructCGALRelated();

	// Store the index of the obstacle vertices, for fast query
	void registerVertices();

	// CGAL and 2D Visibility region
	constrained_delaunay_triangulation::CDT cdt_;
	std::unordered_map<long long, int> cdt_table_;
	int cdt_ver_num_;
	std::vector<std::vector<std::pair<int, int> > > facets_;

	inline int locateAdjacentFacet(std::pair<int, int> prev, std::pair<int, int> next)
	{
		return cdt_table_.find((long long)(prev.first + prev.second * xsize_) + (long long)(next.first + next.second * xsize_) * (long long)(xsize_*ysize_))->second;
	}

};
} // end namespace

