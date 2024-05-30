#pragma once
#define WITH_ROS true
#define WITH_TIME_EVALUATION true
#define WITH_OPENCV_VISUALIZATION false

#include <boost/thread.hpp>

#include <costmap_2d/costmap_2d.h>
#include <ros/ros.h>
#include <costmap_2d/costmap_2d_ros.h>
#include <nav_core/base_global_planner.h>
#include <geometry_msgs/PoseStamped.h>
#include <nav_msgs/GetPlan.h>
#include <angles/angles.h>
#include <base_local_planner/world_model.h>
#include <base_local_planner/costmap_model.h>

namespace raystar 
{
// Check whether the query point is in or on the given polygon
// Return value: <in, on>
std::pair<bool, bool> pnpoly(int nvert, double* vertx, double* verty, 
	double testx, double testy);
std::pair<bool, bool> pnpoly(const std::vector<std::pair<double, double> >& ver, 
	double testx, double testy);

class Candidate
{
public: 
	int Nindex_;
	int Cindex_;
	double Fcost_;

	Candidate(const Candidate& a)
	{
		Nindex_ = a.Nindex_;
		Cindex_ = a.Cindex_;
		Fcost_ = a.Fcost_;
	}
	Candidate(int node_index, int child_index, double cost)
	{
		Nindex_ = node_index;
		Cindex_ = child_index;
		Fcost_ = cost;
	}
};

class Child
{
public:
	// the index of its node
	int Nindex_;

	// the index of itself
	int Cindex_;

	double start_angle_;
	double end_angle_;

	// The "c"
	std::pair<int, int> c_;// The obstacle vertex that form the gap

	// The "o"
	std::pair<double, double> o_;

	int c_obs_index_;
	int c_ver_index_;
	int o_obs_index_;
	int o_ver_index_;

	bool is_a_left_gap_;

	double c_gcost_; // The accumulated distance from the starting location to c
	double c_hcost_; // The heuristic cost

	Child(int nindex, int cindex, int cx, int cy, 
		bool is_left_gap)
	{
		Nindex_ = nindex;
		Cindex_ = cindex;
		c_.first = cx;
		c_.second = cy;
		is_a_left_gap_ = is_left_gap;
	}

	Child(const Child& a)
	{
		Nindex_ = a.Nindex_;
		Cindex_ = a.Cindex_;
		start_angle_ = a.start_angle_;
		end_angle_ = a.end_angle_;
		c_ = a.c_;
		o_ = a.o_;

		c_obs_index_ = a.c_obs_index_;
		c_ver_index_ = a.c_ver_index_;
		o_obs_index_ = a.o_obs_index_;
		o_ver_index_ = a.o_ver_index_;

		is_a_left_gap_ = a.is_a_left_gap_;
		c_gcost_ = a.c_gcost_;
		c_hcost_ = a.c_hcost_;
	}
};

class Polymap;
class Node
{
public:

	int Nindex_;

	// The source point
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

	std::vector< std::pair<double, double> > V_;
	std::vector<std::pair<int, int> > topo_V_;

	std::vector<std::pair<int, int> > local_shortest_path_;
	std::vector<int> path_node_index_;
	
	Node(const Polymap* pMap, int Nindex, double start_x, double start_y, double Gcost, double Hcost, 
		const std::vector<std::pair<double, double> >& visibility_region, 
		const std::vector<std::pair<int, int> >& topo_V);

	Node(const Polymap* pMap, int Nindex, double seed_x, double seed_y, double Gcost, double Hcost,
		int parent_index,
		const std::vector<std::pair<double, double> >& visibility_region,
		const std::vector<std::pair<int, int> >& topo_V);

	void generateChild(const Polymap* pMap);
};

struct PathSolution
{
	std::vector<std::pair<int, int> > path_;
	double path_cost_;
	std::vector<int> path_node_index_;
	PathSolution(const std::vector<std::pair<int, int> >& the_path,
		double the_path_cost, std::vector<int> the_path_node_index)
	{
		path_.assign(the_path.begin(), the_path.end());
		path_cost_ = the_path_cost;
		path_node_index_.assign(the_path_node_index.begin(), the_path_node_index.end());
	}
};

class Raystar: public nav_core::BaseGlobalPlanner
{
public:

    Raystar();
    Raystar(std::string name, costmap_2d::Costmap2D* costmap, std::string frame_id);

    ~Raystar();

    /** overridden classes from interface nav_core::BaseGlobalPlanner **/
    void initialize(std::string name, costmap_2d::Costmap2DROS* costmap_ros);

    void initialize(std::string name, costmap_2d::Costmap2D* costmap, std::string frame_id);

    bool makePlan(const geometry_msgs::PoseStamped& start,
                const geometry_msgs::PoseStamped& goal,
                std::vector<geometry_msgs::PoseStamped>& plan);

	void getScopedVisibilityRegion(Polymap& theMap, 
		Candidate& the_child,
		std::vector<std::pair<double, double> >& visibility_region,
		std::vector<std::pair<int, int> >& topo_V);

	// Storage of computational time, for evaluation
	double map_time_ms_;
	double plan_time_ms_;

	// Storage of the k-SNPP
	std::vector<PathSolution> path_solutions_;

	bool makePlanService(nav_msgs::GetPlan::Request& req, nav_msgs::GetPlan::Response& resp);

protected: 

    costmap_2d::Costmap2D* costmap_;

	// The frame of the map
    std::string frame_id_;

	// Compatiable with existing interface
    ros::Publisher plan_pub_;

	// Publish k-SNPP
	ros::Publisher non_homotopic_pub_;

	// Publish the modified obstacles, for visualization
    ros::Publisher poly_obstacle_pub_;

    bool initialized_;
	
	bool allow_unknown_;

	// Number of required non-homotopic paths
    int K_;

private:
    void clearRobotCell(const geometry_msgs::PoseStamped& global_pose, unsigned int mx, unsigned int my);

    void outlineMap(unsigned char* costarr, int nx, int ny, unsigned char value);

	std::vector<Candidate> Q_; // Open list
	std::vector<Node> N_; // Node list

    boost::mutex mutex_;

    ros::ServiceServer make_plan_srv_;
};
}

