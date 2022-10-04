#include "handle.h"

#include <message_filters/subscriber.h>
#include <message_filters/synchronizer.h>
#include <message_filters/sync_policies/approximate_time.h>

namespace BehaviorDecision{
Handle::Handle(ros::NodeHandle &nodeHandle){
    ROS_INFO("Constructing Handle");
    ROS_INFO("Loading Handle Parameters");

    std::string ego_state_topic_name, rline_topic_name, obstacles_topic_name;
    std::string trajectory_topic_name, viz_trajectory_topic_name;
    if (!nodeHandle.param<std::string>("ego_state_topic_name", ego_state_topic_name,"/ego_state"))
        ROS_WARN_STREAM("Cannot load ego_state_topic_name. Standard value is: " << ego_state_topic_name);
    if (!nodeHandle.param<std::string>("rline_topic_name", rline_topic_name,"/rline"))
        ROS_WARN_STREAM("Cannot load rline_topic_name. Standard value is: " << rline_topic_name);
    if (!nodeHandle.param<std::string>("obstacles_topic_name", obstacles_topic_name,"/obstacles"))
        ROS_WARN_STREAM("Cannot load obstacles_topic_name. Standard value is: " << obstacles_topic_name);
    if (!nodeHandle.param<std::string>("trajectory_topic_name", trajectory_topic_name,"/trajectory/dp_line"))
        ROS_WARN_STREAM("Cannot load trajectory_topic_name. Standard value is: " << trajectory_topic_name);

    if (!nodeHandle.param("node_rate", node_rate, 10))
        ROS_WARN_STREAM("Cannot load node_rate. Standard value is: " << node_rate);
    
    if (!nodeHandle.param("viz", is_viz, true))
        ROS_WARN_STREAM("Cannot load is_viz. Standard value is: " << is_viz);
    if (!nodeHandle.param<std::string>("viz_trajectory_topic_name", viz_trajectory_topic_name,"/viz/trajectory/dp_line"))
        ROS_WARN_STREAM("Cannot load viz_trajectory_topic_name. Standard value is: " << viz_trajectory_topic_name);

    ROS_INFO("Subscribe to Topics");
    message_filters::Subscriber<msgs::Object> ego_state_sub;
    message_filters::Subscriber<msgs::ObjectList> obstacles_sub;
    typedef message_filters::sync_policies::ApproximateTime<msgs::Object, msgs::ObjectList> syncPolicy;
    typedef message_filters::Synchronizer<syncPolicy> Sync;
    boost::shared_ptr<Sync> sync;
    ego_state_sub.subscribe(nodeHandle, ego_state_topic_name,1);
    obstacles_sub.subscribe(nodeHandle, obstacles_topic_name,1);
    sync.reset(new Sync(syncPolicy(10), ego_state_sub, obstacles_sub));   
    sync->registerCallback(boost::bind(&Handle::callback,this, _1, _2));

    ros::Subscriber rline_sub = nodeHandle.subscribe(rline_topic_name, 1, &Handle::rline_callback, this);

    ROS_INFO("Publish to Topics");
    trajectory_pub = nodeHandle.advertise<msgs::Trajectory>(trajectory_topic_name, 1);
    if(is_viz){
        viz_trajectory_pub = nodeHandle.advertise<nav_msgs::Path>(viz_trajectory_topic_name, 1);
    }
    PathP.set_plan_dt(0.1);

}

int Handle::get_NodeRate(){return node_rate;}

void Handle::rline_callback(const msgs::ReferenceLine& rline_data){PathP.set_rline(rline_data);}

void Handle::callback(const msgs::ObjectConstPtr& ego_data, const msgs::ObjectListConstPtr& obstacles_data){
    PathP.path_planning(*ego_data, *obstacles_data);
    VelocityP.save_msgs(PathP.get_msgs());
    VelocityP.velocity_planning();
    }
}

