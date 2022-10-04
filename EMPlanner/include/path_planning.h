#include <ros/ros.h>
#include <msgs/Trajectory.h>
#include <msgs/ReferenceLine.h>
#include <msgs/ObjectList.h>
#include <msgs/Object.h>
#include <map>
#include "velocity_planning.h"

class PathP{
    public:
        PathP();
        void get_param();
        void set_dt(double planning_delta_time);
        void get_referline(const msgs::ReferenceLine rline_data);
        void get_trajactory(const msgs::Trajectory trajectory_data);
        void set_ego_state(const msgs::Object ego_data);
        void set_obstacles(const msgs::ObjectList obstacles_data);

        void calc_startpoint_stitch_trajectory();
        void index2s();
        void dp();

        void get_boundary();
        void qp();
        void path_planning();
        void get_msgs();

        friend class VelocityP;

    private:
        msgs::Object ego_state;
        msgs::ObjectList obstacles;
        msgs::ReferenceLine rline;
        msgs::TrajectoryPoint start_point;
        msgs::Trajectory trajectory;
        msgs::Trajectory pre_trajectory;
        msgs::Trajectory trajectory_dp;
        msgs::Trajectory trajectory_qp;

        msgs::ObjectList static_obstacles;
        msgs::ObjectList dynamic_obstacles;

        std::map<int, int> obs_id2index;

        bool first_run;
        double plan_dt;

        std::vector<double> index2s;
};