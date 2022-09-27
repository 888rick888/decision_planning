#include <ros/ros.h>
#include <msgs/Trajectory.h>
#include <msgs/ReferenceLine.h>
#include <msgs/ObjectList.h>
#include <msgs/Object.h>
#include <map>

class PathP{
    public:
        DP();
        void calc_startpoint_stitch_trajectory();
        void world2frenet_path();
        void calc_dot_infrenet();
        void calc_dot2_infrenet();
        void dp();
        void CalcCost();
        void add_density();
        void path_frenet2cartesian();

        void get_param();
        void set_dt(double planning_delta_time);
        void get_referline(const msgs::ReferenceLine rline_data);
        void get_trajactory(const msgs::Trajectory trajectory_data);
        void set_ego_state(const msgs::Object ego_data);
        void set_obstacles(const msgs::ObjectList obstacles_data);

        void get_boundary()
        void find_near_index()
        void qp()

    private:
        msgs::Object ego_state;
        msgs::ObjectList obstacles;
        msgs::TrajectoryPoint start_point;
        msgs::Trajectory trajectory;
        msgs::Trajectory pre_trajectory;
        msgs::ReferenceLine rline;

        msgs::Trajectory trajectory_qp;

        msgs::ObjectList static_obstacles;
        msgs::ObjectList dynamic_obstacles;

        std::map<int, int> obs_id2index;

        bool first_run;
        double plan_dt;

}