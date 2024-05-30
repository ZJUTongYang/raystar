#include <raystar/raystar.h>
#include <raystar/polymap.h>
#include <chrono>

#include <pluginlib/class_list_macros.h>
#include <visualization_msgs/Marker.h>
#include <visualization_msgs/MarkerArray.h>
#include <angles/angles.h>

PLUGINLIB_EXPORT_CLASS(raystar::Raystar, nav_core::BaseGlobalPlanner)

namespace raystar 
{
std::pair<bool, bool> pnpoly(int nvert, double* vertx, double* verty, double testx, double testy)
{
    std::pair<bool, bool> result; // in, on
    int i, j;
    bool rightray = false;
    bool leftray = false;
    double temp;
    for (i = 0, j = nvert - 1; i < nvert; j = i++) {
        temp = (vertx[j] - vertx[i]) * (testy - verty[i]) / (verty[j] - verty[i]) + vertx[i];
        if (((verty[i] > testy) != (verty[j] > testy)) &&
            (testx < temp))
            rightray = !rightray;
        if (((verty[i] > testy) != (verty[j] > testy)) &&
            (testx > temp))
            leftray = !leftray;
    }
    if (rightray == leftray)
    {
        result.first = rightray;
        result.second = false;
    }
    else
    {
        result.first = false;
        result.second = true;
    }
    return result;
}

std::pair<bool, bool> pnpoly(const std::vector<std::pair<double, double> >& ver, double testx, double testy)
{
    std::pair<bool, bool> result; // in, on
    int i, j;
    bool rightray = false;
    bool leftray = false;
    for (i = 0, j = ver.size() - 1; i < ver.size(); j = i++) {
        if (((ver[i].second > testy) != (ver[j].second > testy)) &&
            (testx < (ver[j].first - ver[i].first) * (testy - ver[i].second) / (ver[j].second - ver[i].second) + ver[i].first))
            rightray = !rightray;
        if (((ver[i].second > testy) != (ver[j].second > testy)) &&
            (testx > (ver[j].first - ver[i].first) * (testy - ver[i].second) / (ver[j].second - ver[i].second) + ver[i].first))
            leftray = !leftray;
    }
    if (rightray == leftray)
    {
        result.first = rightray;
        result.second = false;
    }
    else
    {
        result.first = false;
        result.second = true;
    }
    return result;
}

Node::Node(const Polymap* pMap, int Nindex, double seed_x, double seed_y, double Gcost, double Hcost, 
    const std::vector<std::pair<double, double> >& visibility_region, 
    const std::vector<std::pair<int, int> >& topo_V)
{
    seed_.first = seed_x;
    seed_.second = seed_y;
    Nindex_ = Nindex;
    parent_index_ = -1;
    Gcost_ = Gcost;
    Hcost_ = Hcost;
    local_shortest_path_.emplace_back(std::pair<int, int>(seed_x, seed_y));
    path_node_index_.emplace_back(Nindex_);

    V_.assign(visibility_region.begin(), visibility_region.end());
    topo_V_.assign(topo_V.begin(), topo_V.end());

    start_angle_ = 0;
    end_angle_ = 2 * M_PI;
}

Node::Node(const Polymap* pMap, int Nindex, double seed_x, double seed_y, double Gcost, double Hcost,
    int parent_index, 
    const std::vector<std::pair<double, double> >& visibility_region,
    const std::vector<std::pair<int, int> >& topo_V)
{
    seed_.first = seed_x;
    seed_.second = seed_y;
    Nindex_ = Nindex;
    parent_index_ = parent_index;    
    Gcost_ = Gcost;
    Hcost_ = Hcost;

    V_.assign(visibility_region.begin(), visibility_region.end());
    topo_V_.assign(topo_V.begin(), topo_V.end());
}

void Node::generateChild(const Polymap* pMap)
{
    std::vector<double> theta_list;
    theta_list.resize(V_.size());
    for(unsigned int i = 0; i < V_.size(); ++i)
    {
        theta_list[i] = atan2(V_[i].second - seed_.second, V_[i].first - seed_.first);
        theta_list[i] = start_angle_ + angles::normalize_angle_positive(theta_list[i] - start_angle_);
    }

    std::vector<bool> is_a_gap;
    is_a_gap.resize(V_.size(), false);
    std::vector<bool> is_a_left_gap;
    is_a_left_gap.resize(V_.size(), false);

    std::vector<int> valid_child_indices;

    std::pair<int, int> topo_V_i, topo_V_next;
    for (unsigned int i = 0; i < V_.size()-1; ++i)
    {
        int next = i + 1;
        double angle_diff = angles::normalize_angle(theta_list[next] - theta_list[i]);
        constexpr double threshold_squared = 0.0001 * 0.0001;
        if (angle_diff * angle_diff < threshold_squared)
        {
            topo_V_i = topo_V_[i];
            topo_V_next = topo_V_[next];
            if (!pMap->areConsecutive(topo_V_next, topo_V_i))
            {
                is_a_gap[i] = true;

                double disi = (V_[i].first - seed_.first) * (V_[i].first - seed_.first) + (V_[i].second - seed_.second) * (V_[i].second - seed_.second);
                double disnext = (V_[next].first - seed_.first) * (V_[next].first - seed_.first) + (V_[next].second - seed_.second) * (V_[next].second - seed_.second);
                if (disi > disnext)
                    is_a_left_gap[i] = true;
                valid_child_indices.emplace_back(i);
            }
        }
    }

    // For each gap, we construct a child structure
    for(auto iter = valid_child_indices.begin(); iter != valid_child_indices.end(); ++iter)
    {
        int i = *iter;
        int next = (i + 1) % V_.size();
        topo_V_i = topo_V_[i];
        topo_V_next = topo_V_[next];
        if(is_a_left_gap[i])
        {
            C_.emplace_back(Child(Nindex_, -1, (int)(V_[next].first), (int)(V_[next].second), is_a_left_gap[i]));
            
            C_.back().start_angle_ = angles::normalize_angle(theta_list[next]);
            std::pair<double, double> next_obs = pMap->getNextObs(topo_V_next);
            double contour_angle_from_next = atan2(next_obs.second - V_[next].second, 
                next_obs.first - V_[next].first);
            C_.back().end_angle_ = C_.back().start_angle_ + angles::normalize_angle_positive(contour_angle_from_next - C_.back().start_angle_);

            C_.back().c_obs_index_ = topo_V_next.first;
            C_.back().c_ver_index_ = topo_V_next.second;
            C_.back().o_obs_index_ = topo_V_i.first;
            C_.back().o_ver_index_ = topo_V_i.second;
            C_.back().o_.first = V_[i].first;
            C_.back().o_.second = V_[i].second;
        }
        else
        {
            C_.emplace_back(Child(Nindex_, -1, (int)(V_[i].first), (int)(V_[i].second), is_a_left_gap[i]));

            std::pair<double, double> prev_obs = pMap->getPrevObs(topo_V_i);
            double contour_angle_from_prev = atan2(prev_obs.second - V_[i].second, 
                prev_obs.first - V_[i].first);
            C_.back().start_angle_ = contour_angle_from_prev;
            C_.back().end_angle_ = contour_angle_from_prev + angles::normalize_angle_positive(theta_list[i] - contour_angle_from_prev);

            C_.back().c_obs_index_ = topo_V_i.first;
            C_.back().c_ver_index_ = topo_V_i.second;
            C_.back().o_obs_index_ = topo_V_next.first;
            C_.back().o_ver_index_ = topo_V_next.second;
            C_.back().o_.first = V_[next].first;
            C_.back().o_.second = V_[next].second;
        }
    }

    // After sorting the children based on gap orientation, 
    for (unsigned int i = 0; i < C_.size(); ++i)
    {
        C_[i].Cindex_ = i;
    }

    for (auto iter = C_.begin(); iter != C_.end(); ++iter)
    {
        iter->c_gcost_ = Gcost_ + sqrt((seed_.first - iter->c_.first) * (seed_.first - iter->c_.first) + (seed_.second - iter->c_.second) * (seed_.second - iter->c_.second));
    }

}

void Raystar::outlineMap(unsigned char* costarr, int nx, int ny, unsigned char value) 
{
    unsigned char* pc = costarr;
    for (int i = 0; i < nx; i++)
        *pc++ = value;
    pc = costarr + (ny - 1) * nx;
    for (int i = 0; i < nx; i++)
        *pc++ = value;
    pc = costarr;
    for (int i = 0; i < ny; i++, pc += nx)
        *pc = value;
    pc = costarr + nx - 1;
    for (int i = 0; i < ny; i++, pc += nx)
        *pc = value;
}

Raystar::Raystar():
    costmap_(NULL), allow_unknown_(true), initialized_(false)
{

}

Raystar::Raystar(std::string name, costmap_2d::Costmap2D* costmap, std::string frame_id):
    Raystar()
{
   initialize(name, costmap, frame_id);
}

Raystar::~Raystar()
{

}

void Raystar::initialize(std::string name, costmap_2d::Costmap2DROS* costmap_ros){
    initialize(name, costmap_ros->getCostmap(), costmap_ros->getGlobalFrameID());
}

void Raystar::initialize(std::string name, costmap_2d::Costmap2D* costmap, std::string frame_id)
{
    if (!initialized_) {
        ros::NodeHandle private_nh("~/" + name);
        costmap_ = costmap;
        frame_id_ = frame_id;

        plan_pub_ = private_nh.advertise<nav_msgs::Path>("plan", 1);

        non_homotopic_pub_ = private_nh.advertise<visualization_msgs::MarkerArray>("non_homotopic_paths", 1);

        poly_obstacle_pub_ = private_nh.advertise<visualization_msgs::MarkerArray>("poly_obstacles", 1);

        private_nh.param("allow_unknown", allow_unknown_, true);

        private_nh.param("num_of_holonomic_paths", K_, 5);

        make_plan_srv_ = private_nh.advertiseService("make_plan", &Raystar::makePlanService, this);

        initialized_ = true;
    } else
        ROS_WARN("This planner has already been initialized, you can't call it twice, doing nothing");
}

void Raystar::clearRobotCell(const geometry_msgs::PoseStamped& global_pose, unsigned int mx, unsigned int my) {
    if (!initialized_) {
        ROS_ERROR(
                "This planner has not been initialized yet, but it is being used, please call initialize() before use");
        return;
    }

    //set the associated costs in the cost map to be free
    costmap_->setCost(mx, my, costmap_2d::FREE_SPACE);
}

bool Raystar::makePlanService(nav_msgs::GetPlan::Request& req, nav_msgs::GetPlan::Response& resp) {
    // For consistency, we only use 1-SNPP for services
    int temp_K = K_;
    K_ = 1;
    makePlan(req.start, req.goal, resp.plan.poses);
    K_ = temp_K;
    
    resp.plan.header.stamp = ros::Time::now();
    resp.plan.header.frame_id = frame_id_;

    return true;
}

void Raystar::getScopedVisibilityRegion(
    Polymap& theMap, 
    Candidate& the_child, 
    std::vector<std::pair<double, double> >& visibility_region, 
    std::vector<std::pair<int, int> >& topo_V)
{
    visibility_region.clear();
    topo_V.clear();

    int parent_index = the_child.Nindex_;
    int child_index = the_child.Cindex_;
    std::pair<int, int> new_source_point = N_[parent_index].C_[child_index].c_;
    std::vector<std::pair<double, double> > fullV;
    std::vector<std::pair<int, int> > full_topoV;
    theMap.getVisibilityRegion(new_source_point.first, new_source_point.second, fullV, full_topoV);

    std::pair<double, double> start_obs;
    std::pair<double, double> end_obs;
    std::pair<int, int> start_obs_topo;
    std::pair<int, int> end_obs_topo;

    double start_angle;
    double end_angle;

    if (N_[parent_index].C_[child_index].is_a_left_gap_)
    {
        start_obs = N_[parent_index].C_[child_index].o_;
        end_obs = theMap.getNextObs(std::pair<int, int>( N_[parent_index].C_[child_index].c_obs_index_, N_[parent_index].C_[child_index].c_ver_index_ ));
        start_obs_topo.first = N_[parent_index].C_[child_index].o_obs_index_;
        start_obs_topo.second = N_[parent_index].C_[child_index].o_ver_index_;
        end_obs_topo = theMap.locateVertex(end_obs.first, end_obs.second);
        start_angle = atan2(start_obs.second - new_source_point.second, start_obs.first - new_source_point.first);
        end_angle = atan2(end_obs.second - new_source_point.second, end_obs.first - new_source_point.first);
    }
    else
    {
        start_obs = theMap.getPrevObs(std::pair<int, int>(N_[parent_index].C_[child_index].c_obs_index_, N_[parent_index].C_[child_index].c_ver_index_));
        end_obs = N_[parent_index].C_[child_index].o_;
        start_obs_topo = theMap.locateVertex(start_obs.first, start_obs.second);
        end_obs_topo.first = N_[parent_index].C_[child_index].o_obs_index_;
        end_obs_topo.second = N_[parent_index].C_[child_index].o_ver_index_;
        start_angle = atan2(start_obs.second - new_source_point.second, start_obs.first - new_source_point.first);
        end_angle = atan2(end_obs.second - new_source_point.second, end_obs.first - new_source_point.first);
    }
    end_angle = start_angle + angles::normalize_angle_positive(end_angle - start_angle);

    double theta;
    for (unsigned int i = 0; i < fullV.size(); ++i)
    {
        theta = atan2(fullV[i].second - new_source_point.second, fullV[i].first - new_source_point.first);
        theta = start_angle + angles::normalize_angle(theta - start_angle);
        if (theta >= start_angle - 0.0000001 && theta <= end_angle + 0.0000001)
        {
            visibility_region.emplace_back(fullV[i]);
            topo_V.emplace_back(full_topoV[i]);
        }
    }

    // We add the "limit"
    int loc = std::find_if(visibility_region.begin(), visibility_region.end(),
        [start_obs](const auto& a) {return a.first == start_obs.first && a.second == start_obs.second; }) - visibility_region.begin();

    if(loc == visibility_region.size())
    {
        visibility_region.insert(visibility_region.begin(), start_obs);
        topo_V.insert(topo_V.begin(), theMap.locateVertex(start_obs.first, start_obs.second));
    }
    else
    {
        visibility_region.erase(visibility_region.begin(), visibility_region.begin() + loc);
        topo_V.erase(topo_V.begin(), topo_V.begin() + loc);
    }

    loc = std::find_if(visibility_region.begin(), visibility_region.end(),
        [end_obs](const auto& a) {return a.first == end_obs.first && a.second == end_obs.second; }) - visibility_region.begin();

    if(loc == visibility_region.size())
    {
        visibility_region.insert(visibility_region.end(), end_obs);
        topo_V.insert(topo_V.end(), theMap.locateVertex(end_obs.first, end_obs.second));
    }
    else
    {
        visibility_region.erase(visibility_region.begin() + loc + 1, visibility_region.end());
        topo_V.erase(topo_V.begin() + loc + 1, topo_V.end());
    }
}

bool Raystar::makePlan(const geometry_msgs::PoseStamped& start, const geometry_msgs::PoseStamped& goal, std::vector<geometry_msgs::PoseStamped>& plan) 
{

    boost::mutex::scoped_lock lock(mutex_);
    if (!initialized_) {
        ROS_ERROR(
            "This planner has not been initialized yet, but it is being used, please call initialize() before use");
        return false;
    }

    //clear the plan, just in case
    plan.clear();
    path_solutions_.clear();

    ros::NodeHandle n;
    std::string global_frame = frame_id_;

    //until tf can handle transforming things that are way in the past... we'll require the goal to be in our global frame
    if (goal.header.frame_id != global_frame) {
        ROS_ERROR(
            "The goal pose passed to this planner must be in the %s frame.  It is instead in the %s frame.", global_frame.c_str(), goal.header.frame_id.c_str());
        return false;
    }

    if (start.header.frame_id != global_frame) {
        ROS_ERROR(
            "The start pose passed to this planner must be in the %s frame.  It is instead in the %s frame.", global_frame.c_str(), start.header.frame_id.c_str());
        return false;
    }

    double wx = start.pose.position.x;
    double wy = start.pose.position.y;

    unsigned int start_x_i, start_y_i, goal_x_i, goal_y_i;
    double start_x, start_y, goal_x, goal_y;

    if (!costmap_->worldToMap(wx, wy, start_x_i, start_y_i)) {
        ROS_WARN_THROTTLE(1.0,
            "The robot's start position is off the global costmap. Planning will always fail, are you sure the robot has been properly localized?");
        return false;
    }
    start_x = start_x_i;
    start_y = start_y_i;

    wx = goal.pose.position.x;
    wy = goal.pose.position.y;

    if (!costmap_->worldToMap(wx, wy, goal_x_i, goal_y_i)) {
        ROS_WARN_THROTTLE(1.0,
            "The goal sent to the global planner is off the global costmap. Planning will always fail to this goal.");
        return false;
    }
    goal_x = goal_x_i;
    goal_y = goal_y_i;

    //clear the starting cell within the costmap because we know it can't be an obstacle
    clearRobotCell(start, start_x_i, start_y_i);

    int nx = costmap_->getSizeInCellsX(), ny = costmap_->getSizeInCellsY();

    outlineMap(costmap_->getCharMap(), nx, ny, costmap_2d::LETHAL_OBSTACLE);

    // We need to obtain a topological representation
#if WITH_TIME_EVALUATION
    auto map_start_time = std::chrono::high_resolution_clock::now();
#endif
    Polymap theMap(costmap_, start_x, start_y, goal_x, goal_y);
#if WITH_TIME_EVALUATION
    auto map_end_time = std::chrono::high_resolution_clock::now();
    auto map_duration = std::chrono::duration_cast<std::chrono::microseconds>(map_end_time - map_start_time);
    map_time_ms_ = map_duration.count() / 1000.0;
#endif

    if(!theMap.solution_exist_)
        return false;

    if (poly_obstacle_pub_.getNumSubscribers() > 0)
    {
        static visualization_msgs::MarkerArray array;
        // We first clear dated markers
        if(array.markers.size() != 0)
        {
            for(auto iter = array.markers.begin(); iter != array.markers.end(); ++iter)
            {
                iter->action = visualization_msgs::Marker::DELETE;
            }
            poly_obstacle_pub_.publish(array);
        }

        // We visualize polygonal obstacles
        static visualization_msgs::Marker marker;
        marker.header.frame_id = global_frame;
        marker.header.stamp = ros::Time::now();
        marker.ns = "polygons";
        marker.id = 0;
        marker.type = visualization_msgs::Marker::LINE_STRIP;
        marker.action = visualization_msgs::Marker::ADD;
        marker.pose.position.x = 0.0;
        marker.pose.position.y = 0.0;
        marker.pose.position.z = 0.0;
        marker.pose.orientation.x = 0.0;
        marker.pose.orientation.y = 0.0;
        marker.pose.orientation.z = 0.0;
        marker.pose.orientation.w = 1.0;
        marker.scale.x = 0.1; // Line width
        marker.color.r = 1.0;
        marker.color.g = 0.0;
        marker.color.b = 0.0;
        marker.color.a = 1.0;
        marker.lifetime = ros::Duration(0);

        for (auto iter = theMap.obs_.begin(); iter != theMap.obs_.end(); ++iter)
        {
            marker.id++;
            marker.points.clear();
            marker.colors.clear();
            geometry_msgs::Point temp;
            std_msgs::ColorRGBA color_temp;
            color_temp.r = 1.0;
            double wx, wy;
            for (auto iter2 = iter->detailed_ordered_vertices_.begin(); iter2 != iter->detailed_ordered_vertices_.end(); ++iter2)
            {
                auto next = std::next(iter2);
                if (next == iter->detailed_ordered_vertices_.end())
                    next = iter->detailed_ordered_vertices_.begin();

                costmap_->mapToWorld(iter2->first, iter2->second, wx, wy);
                temp.x = wx;
                temp.y = wy;
                temp.z = 0;
                marker.points.emplace_back(temp);
                marker.colors.emplace_back(color_temp);
                costmap_->mapToWorld(next->first, next->second, wx, wy);
                temp.x = wx;
                temp.y = wy;
                temp.z = 0;
                marker.points.emplace_back(temp);
                marker.colors.emplace_back(color_temp);

            }
            array.markers.emplace_back(marker);
        }
        poly_obstacle_pub_.publish(array);
    }

#if WITH_TIME_EVALUATION
    auto planner_start_time = std::chrono::high_resolution_clock::now();
#endif
    // We initialise the priority queue
    auto comp = [](const Candidate& a, const Candidate& b) {return a.Fcost_ > b.Fcost_; };

    Q_.clear();
    N_.clear();
    Q_.emplace_back(Candidate(-1, -1, sqrt((start_x - goal_x) * (start_x - goal_x) + (start_y - goal_y) * (start_y - goal_y))));

    while (!Q_.empty())
    {
        Candidate best_candidate = Q_[0];
        std::pop_heap(Q_.begin(), Q_.end(), comp);
        Q_.pop_back();

        int parent_index = best_candidate.Nindex_;
        int child_index = best_candidate.Cindex_;

        if (parent_index == -1)
        {
            // For the first node
            std::vector<std::pair<double, double> > Vtemp;
            std::vector<std::pair<int, int> > topo_Vtemp;
            theMap.getVisibilityRegion(start_x, start_y, Vtemp, topo_Vtemp);
            N_.emplace_back(Node(&theMap, 0, start_x, start_y, 0.0, best_candidate.Fcost_, Vtemp, topo_Vtemp));
            N_.back().generateChild(&theMap);
        }
        else
        {
            // For other non-root nodes
            std::pair<double, double> new_source_point = N_[parent_index].C_[child_index].c_;
            int new_node_index = N_.size();
            std::vector<std::pair<double, double> > Vtemp;
            std::vector<std::pair<int, int> > topo_Vtemp;
            getScopedVisibilityRegion(theMap, best_candidate, Vtemp, topo_Vtemp);

            N_.emplace_back(Node(&theMap, new_node_index, new_source_point.first, new_source_point.second,
                N_[parent_index].C_[child_index].c_gcost_, N_[parent_index].C_[child_index].c_hcost_,
                parent_index, Vtemp, topo_Vtemp));

            N_.back().local_shortest_path_.assign(N_[parent_index].local_shortest_path_.begin(), N_[parent_index].local_shortest_path_.end());
            N_.back().local_shortest_path_.emplace_back(new_source_point);
            N_.back().path_node_index_.emplace_back(new_node_index);
            N_.back().start_angle_ = N_[parent_index].C_[child_index].start_angle_;
            N_.back().end_angle_ = N_[parent_index].C_[child_index].end_angle_;

            N_.back().generateChild(&theMap);
        }

        for (auto iter = N_.back().C_.begin(); iter != N_.back().C_.end(); ++iter)
        {
            iter->c_hcost_ = sqrt((iter->c_.first - goal_x) * (iter->c_.first - goal_x) + (iter->c_.second - goal_y) * (iter->c_.second - goal_y));

            Q_.emplace_back(Candidate(iter->Nindex_, iter->Cindex_, iter->c_gcost_ + iter->c_hcost_));
            std::push_heap(Q_.begin(), Q_.end(), comp);
        }

        // We check whether the goal is within the visibility region
        std::pair<bool, bool> b = pnpoly(N_.back().V_, goal_x, goal_y);
        if (b.first || b.second)
        {
            std::vector<std::pair<int, int> > locally_shortest_path(N_.back().local_shortest_path_);
            locally_shortest_path.emplace_back(std::pair<int, int>(goal_x, goal_y));
            double new_path_length = N_.back().Gcost_ + sqrt((N_.back().seed_.first - goal_x) * (N_.back().seed_.first - goal_x)
                + (N_.back().seed_.second - goal_y) * (N_.back().seed_.second - goal_y));

            path_solutions_.emplace_back(PathSolution(locally_shortest_path, new_path_length, N_.back().path_node_index_));
            if (path_solutions_.size() >= K_)
            {
                break;
            }
        }
    }
#if WITH_TIME_EVALUATION
    auto planner_end_time = std::chrono::high_resolution_clock::now();
    auto planner_duration = std::chrono::duration_cast<std::chrono::microseconds>(planner_end_time - planner_start_time);
    plan_time_ms_ = planner_duration.count() / 1000.0;
#endif


    if (non_homotopic_pub_.getNumSubscribers() > 0)
    {
        static visualization_msgs::MarkerArray array;
        // We first clear dated markers
        if(array.markers.size() != 0)
        {
            for(auto iter = array.markers.begin(); iter != array.markers.end(); ++iter)
            {
                iter->action = visualization_msgs::Marker::DELETE;
            }
            non_homotopic_pub_.publish(array);
        }

        // We create a proper colour table for visualization
        int num_for_division = ceil( sqrt(path_solutions_.size()) );
        int step = floor(100 / (num_for_division+1));
        std::vector<std::vector<int> > colortable;
        std::vector<int> colortemp(3);
        colortable.resize((num_for_division+1)*(num_for_division+1)*(num_for_division+1), colortemp);

        for(unsigned int i = 0; i < num_for_division; ++i)
        {
            for(unsigned int j = 0; j < num_for_division; ++j)
            {
                for(unsigned int k = 0; k < num_for_division; ++k)
                {
                    colortable[i*(num_for_division+1)*(num_for_division+1)+j*(num_for_division+1)+k] = {100+i*step, 100+j*step, 100+k*step};
                }
            }
        }

        visualization_msgs::Marker marker;
        marker.header.frame_id = global_frame;
        marker.header.stamp = ros::Time::now();
        marker.ns = "non_homotopic_paths";
        marker.id = 0;
        marker.type = visualization_msgs::Marker::LINE_STRIP;
        marker.action = visualization_msgs::Marker::ADD;
        marker.pose.position.x = 0.0;
        marker.pose.position.y = 0.0;
        marker.pose.position.z = 0.0;
        marker.pose.orientation.x = 0.0;
        marker.pose.orientation.y = 0.0;
        marker.pose.orientation.z = 0.0;
        marker.pose.orientation.w = 1.0;
        marker.scale.x = 0.1; // Line width
        marker.color.r = 0.0;
        marker.color.g = 0.0;
        marker.color.b = 0.0;
        marker.color.a = 1.0;
        marker.lifetime = ros::Duration(0);

        for(auto iter = path_solutions_.begin(); iter != path_solutions_.end(); ++iter)
        {

            marker.id++;
            marker.points.clear();
            marker.colors.clear();
            geometry_msgs::Point temp;
            std_msgs::ColorRGBA color_temp;
            color_temp.r = colortable[iter - path_solutions_.begin()][0] / 255.0;
            color_temp.g = colortable[iter - path_solutions_.begin()][1] / 255.0;
            color_temp.b = colortable[iter - path_solutions_.begin()][2] / 255.0;
            double wx, wy;
            for (auto iter2 = iter->path_.begin(); iter2 != iter->path_.end()-1; ++iter2)
            {
                auto next = std::next(iter2);

                costmap_->mapToWorld(iter2->first, iter2->second, wx, wy);
                temp.x = wx;
                temp.y = wy;
                temp.z = 0;
                marker.points.emplace_back(temp);
                marker.colors.emplace_back(color_temp);
                costmap_->mapToWorld(next->first, next->second, wx, wy);
                temp.x = wx;
                temp.y = wy;
                temp.z = 0;
                marker.points.emplace_back(temp);
                marker.colors.emplace_back(color_temp);
            }
            array.markers.emplace_back(marker);
        }
        non_homotopic_pub_.publish(array);
    }

    if(!path_solutions_.empty())
    {
        // We create the path in the world frame
        std::vector<std::pair<double, double> > path_double;
        path_double.reserve(ceil(path_solutions_[0].path_cost_));
        for(unsigned int i = 0; i < path_solutions_[0].path_.size()-1; ++i)
        {
            std::pair<double, double> prev = path_solutions_[0].path_[i];
            std::pair<double, double> next = path_solutions_[0].path_[i+1];
            double dis = sqrt( (prev.first-next.first)*(prev.first-next.first) + (prev.second - next.second)*(prev.second - next.second) );
            int count = ceil(dis);
            for(unsigned int j = 0; j < count; ++j)
            {
                path_double.emplace_back(std::pair<double, double>( 
                    (prev.first*(count-j)+next.first*j)/count, (prev.second*(count-j)+next.second*j)/count
                 ));
            }
        }
        path_double.emplace_back(path_solutions_[0].path_.back());

        geometry_msgs::PoseStamped temp;
        temp.header.stamp = ros::Time::now();
        temp.header.frame_id = frame_id_;
        temp.pose.orientation.w = 1;
        for(auto iter = path_double.begin(); iter != path_double.end(); ++iter)
        {
            double world_x, world_y;
            costmap_->mapToWorld(iter->first, iter->second, world_x, world_y);
            temp.pose.position.x = world_x;
            temp.pose.position.y = world_y;
            plan.emplace_back(temp);
        }
    }

#if WITH_TIME_EVALUATION
    std::cout << "Time for transforming gridmap to polygonal map: " << map_time_ms_ << "ms, time for planning: " << plan_time_ms_ << "ms. " << std::endl;
#endif

    return !plan.empty();
}

};