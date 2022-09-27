#pragma once

#include <ros/ros.h>
#include <msgs/ObjectList.h>
#include <msgs/ReferenceLine.h>
#include <msgs/Trajectory.h>
#include <nav_msgs/Path.h>

#include <dp.h>
namespace BehaviorDecision{
class Handle{
    public:
        Handle(ros::NodeHandle &nodeHandle);
        int get_NodeRate();
        void callback(const msgs::ObjectConstPtr& ego_data, const msgs::ObjectListConstPtr& obstacles_data);
        void rline_callback(const msgs::ReferenceLine& rline);

    private:

        ros::Publisher trajectory_pub;
        ros::Publisher viz_trajectory_pub;
        int node_rate;
        bool is_viz;

        DP dp;
};
}