#include "handle.h"

int main(int argc,char** argv)
{
    ros::init(argc, argv, "behavior_decision");
    ros::NodeHandle nodeHandle("~");
    BehaviorDecision::Handle handle(nodeHandle);
    ros::Rate loop_rate(handle.get_NodeRate());
    while(ros::ok())
    {
        ros::spinOnce();   // Keeps node alive basically
        loop_rate.sleep(); // Sleep for loop_rate
    }
    return 0;
}
