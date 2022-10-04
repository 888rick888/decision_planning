#include <ros/ros.h>
#include <msgs/Trajectory.h>
#include <msgs/ReferenceLine.h>
#include <msgs/ObjectList.h>
#include <msgs/Object.h>
#include <map>
#include <path_planning.h>

class VelocityP{
    public:
        VelocityP();
        void get_param();
        void set_dt(double planning_delta_time);

        void trajectory_s2xy();
        void programming_init();
        void generate_st_graph();
        void speed_dp();
        double CalcDPCost(int row_start, int col_start, int row_end, int col_end, auto s_list, auto t_list, auto dp_st_s_dot);
        void Convex_space();
        void qp_velocity();
        void velocity_planning();
        void save_msgs(msgs::TrajectoryPoint start, msgs::Object ego, msgs::ObjectList dynamic, msgs::Trajectory trajectory);

    private:
        msgs::Object ego_state;
        msgs::TrajectoryPoint start_point;
        msgs::Trajectory trajectory_qp;
        msgs::Trajectory trajectory_init;

        msgs::Velocity_planning_st obs_st;
        msgs::ObjectList dynamic_obstacles;

        std::map<int, int> obs_id2index;

        bool first_run;
        double plan_dt;

        std::vector<double> path_index2s;
};