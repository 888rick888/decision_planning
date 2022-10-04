#include "velocity_planning.h"
#include <map>
#include <cmath>
#include <eigen3/Eigen/Dense>
#include <string>
#include <ros/ros.h>
#include <vector>
#include "OsqpEigen/OsqpEigen.h"
#include <msgs/ReferencePoint.h>
#include <msgs/TrajectoryPoint.h>
#include <msgs/Object.h>
#include <msgs/ObjectList.h>
#include <msgs/Velocity_planning_st.h>
#include <tools/tools.hpp>
#include <tools/frenet.hpp>
#include <tools/velocityp_tools.hpp>

double reference_speed;
int w_cost_ref_speed, w_cost_accel, w_cost_obs;
int w_cost_s_dot2, w_cost_v_ref, w_cost_cost_jerk;

VelocityP::VelocityP(){
    first_run = true;
}

void VelocityP::get_param(){
    ros::param::get("/reference_speed", reference_speed);
    ros::param::get("/w_cost_ref_speed", w_cost_ref_speed);
    ros::param::get("/w_cost_accel", w_cost_accel);
    ros::param::get("/w_cost_obs", w_cost_obs);
    ros::param::get("/w_cost_s_dot2", w_cost_s_dot2);
    ros::param::get("/w_cost_v_ref", w_cost_v_ref);
    ros::param::get("/w_cost_cost_jerk", w_cost_cost_jerk);
}

void VelocityP::get_referline(const msgs::ReferenceLine rline_data){
    rline = rline_data;
}

void VelocityP::set_dt(double planning_delta_time){
    plan_dt = planning_delta_time;
}

void velocity_planning(){
    get_param();
    trajectory_s2xy();
    programming_init();
    generate_st_graph();
    speed_dp();
    Convex_space();
    qp_velocity();
}

void trajectory_s2xy(){     //////////trajectory_init, 为路径规划的转换成世界坐标系的最终结果
    int sum = 0;

    int i = 0;
    for (; i<trajectory_init.points.size(); i++){
        sum += std::sqrt(std::pow(trajectory_init.points[i].pose.x - trajectory_init.points[i-1].pose.x, 2) + std::pow(trajectory_init.points[i].pose.y - trajectory_init.points[i-1].pose.y));
        path_index2s[i] = sum;
    }
    double path_s_end = path_index2s[-1];
}

void programming_init(){        /////////////////start_point
    start_point.frenet_s[1] = std::cos(start_point.pose.theta) * start_point.velocity.x + std::sin(start_point.pose.theta) * start_point.velocity.y;
    start_point.frenet_s[2] = std::cos(start_point.pose.theta) * start_point.accel.x + std::sin(start_point.pose.theta) * start_point.accel.y;
}

void generate_st_graph(){        /////////////////dynamic_obstacles
    for (int i=0; i<dynamic_obstacles.objects.size(); i++){
        if (std::abs(dynamic_obstacles.objects[i].frenet_l[1]) < 0.3){
            if (std::abs(dynamic_obstacles.objects[i].frenet_l[0]) > 2){
                continue;
            }
            else{continue;  //ODO需要做虚拟障碍物，感知模块加入跟踪逻辑，判断两帧之间的障碍物是否是同一个，给到速度规划模块后，先给出虚拟障碍物觉，下一帧拿到虚拟障碍物的信息，规划绕过去的路径和速度
            }               //本算法欠缺的：障碍物结构体，包含了坐标、速度、还要包含决策标记（是否为虚拟障碍物，左绕还是右绕，避让还是超车
        }
        double t_zero = -dynamic_obstacles.objects[i].frenet_l[0] / dynamic_obstacles.objects[i].frenet_l[1];
        double t_boundary1 = 2 / dynamic_obstacles.objects[i].frenet_l[1] + t_zero;
        double t_boundary2 = 2 / dynamic_obstacles.objects[i].frenet_l[1] + t_zero;
        double t_max = std::max(t_boundary1, t_boundary2);
        double t_min = std::min(t_boundary1, t_boundary2);

        if (t_max < 1 || t_min > 8){
            continue;
        }
        if (t_min < 0 && t_max > 0){    //感知看到时，障碍物已经在±2的内部了
            obs_st.s_in[i] = dynamic_obstacles.objects[i].frenet_s[0];
            obs_st.s_out[i] = dynamic_obstacles.objects[i].frenet_s[0] + dynamic_obstacles.objects[i].frenet_s[1] * t_max;
            obs_st.t_in[i] = 0;
            obs_st.t_out[i] = t_max;
        }
        else{       //正常障碍物
            obs_st.s_in[i] = dynamic_obstacles.objects[i].frenet_s[0] + dynamic_obstacles.objects[i].frenet_s[1] * t_min;
            obs_st.s_out[i] = dynamic_obstacles.objects[i].frenet_s[0] + dynamic_obstacles.objects[i].frenet_s[1] * t_max;
            obs_st.t_in[i] = t_min;
            obs_st.t_out[i] = t_max;
        }
    }
}


std::vector<int> dp_speed_s(16, -1);
std::vector<int> dp_speed_t(16, -1);
void speed_dp(){
    Eigen::MatrixXd <double, 40, 16> dp_st_cost = Eigen::Ones();
    Eigen::MatrixXd <double, 40, 16> dp_st_node = Eigen::Ones();
    Eigen::MatrixXd <double, 40, 16> dp_st_s_dot= Eigen::Zeros();   //表示从起点到i，j点的最优路径的末速度
    // dp_st_cost.setOnes();
    dp_st_cost = dp_st_cost * 99999;

    std::vector<double> s_list(40), t_list(16);
    clac_arithmetic(s_list, 0, 4.5, 0.5, 1);
    clac_arithmetic(s_list, 5.5, 14.5, 1, 2);
    clac_arithmetic(s_list, 16, 29.5, 1.5, 3);
    clac_arithmetic(s_list, 32, 54.5, 2.5, 4);
    clac_arithmetic(t_list, 0.5, 8, 0.5, 1);

    for (int i=0; i<s_list.size(); i++){
        dp_st_cost(i, 0) = CalcDPCost(0, 0, i, 0, s_list, t_list, dp_st_s_dot);
        auto s_end, t_end = CalcSTCoordinate(i, 0, s_list, t_list);
        dp_st_s_dot(i, 0) = s_end / t_end;
    }

    for (int i=1, i<t_list.size(); i++){
        for (int j=0; j<s_list.size(); j++){
            int cur_row = j;
            int cur_col = i;
            for (int k=0; K<s_list.size(); k++){
                int pre_row = k;
                int pre_col = i - 1;
                double cost_temp = CalcDPCost(pre_row, pre_col, cur_row, cur_col, s_list, t_list, dp_st_s_dot);
                if (cost_temp + dp_st_cost(pre_row, pre_col) < dp_st_cost(cur_row, cur_col)){
                    dp_st_cost(cur_row, cur_col) = cost_temp + dp_st_cost(pre_row, pre_col);
                    auto s_start, t_start = CalcSTCoordinate(pre_row, pre_col, s_list, t_list);
                    auto s_end, t_end = CalcSTCoordinate(cur_row, cur_col, s_list, t_list);
                    dp_st_s_dot(cur_row, cur_col) = (s_end - s_start) / (t_end - t_start);
                    dp_st_node(cur_row, cur_col) = pre_row; 
                }
            }
        }
    }
    int min_cost = 99999, min_row = 99999, min_col = 99999;
    for (int i=0; i<s_list.size(); i++){
        if (dp_st_cost=(i, -1) < min_cost){
            min_cost = dp_st_cost(i, -1);
            min_row = i;
            min_col = t_list.size();
        }
    }
    for (int i=0; i<t_list.size(); i++){
        if (dp_st_cost(0, j) <= min_cost){
            min_cost = dp_st_cost;
            min_row = 0;
            min_col = i;
        }
    }

    dp_speed_s[min_col], dp_speed_t[min_col] = CalcSTCoordinate(min_row, min_col, s_list, t_list);
    while (min_col != 1){
        pre_row = dp_st_node(min_row, min_col);
        pre_col = min_col - 1;
        dp_speed_s[pre_col], dp_speed_t[pre_col] = CalcSTCoordinate(pre_row, pre_col, s_list, t_list);
        min_row = pre_row;
        min_col = pre_col;
    }
}

double CalcDPCost(int row_start, int col_start, int row_end, int col_end, auto s_list, auto t_list, auto dp_st_s_dot){
    auto s_end, t_end = CalcSTCoordinate(row_end, col_end, s_list, t_list);
    if (row_start == 0){
        int s_start = 0;
        int t_start = 0;
        auto s_dot_start = start_point.frenet_s[1];
    }
    else{
        auto s_start, t_start = CalcSTCoordinate(row_start, col_start, s_list, t_list);
        auto s_dot_start = dp_st_s_dot(row_start, col_start);
    }

    double cur_s_dot = (s_end - s_start)/ (t_end - t_start);
    double cur_s_dot2 = (cur_s_dot - s_dot_start) / (t_end - t_start);
    double cost_ref_speed = w_cost_ref_speed * std::pow((cur_s_dot - reference_speed), 2);
    if (cur_s_dot2 < 4 && cur_s_dot2 > -6){
        cost_accel = w_cost_accel * std::pow(cur_s_dot2, 2);
    }
    else{cost_accel = 9999 * w_cost_accel * std::pow(cur_s_dot2, 2)}

    double cost_obs = CalcObsCost(s_start, t_start, s_end, t_end, obs_st, w_cost_obs);
    cost = cost_obs + cost_accel + cost_ref_speed;
    return cost
}


std::vector s_lb(dp_speed_s.size(), -99999);
std::vector s_ub(dp_speed_s.size(), 99999);
std::vector s_dot_lb(dp_speed_s.size(), -99999);
std::vector s_dot_ub(dp_speed_s.size(), 99999);
void Convex_space(){
    int path_index2s_end_index = path_index2s.size();
    int dp_speed_end_index = dp_speed_s.size();

    // for (int k=1; k<path_index2s.size(); k++){
    //     if (path_index2s[k] == 0 && path_index2s[k-1] != 0){
    //         path_index2s_end_index = k - 1;
    //         break;
    //     }
    //     path_index2s_end_index = k;
    // }

    for (int k=0; k<dp_speed_s.size(); k++){
        if (dp_speed_s[k] < 0){
            dp_speed_end_index = k - 1;
            break;
        }
    }

    for (int i=0; dp_speed_s[i]>0; i++){    //施加车辆动力学约束
        if (dp_speed_s[i] < 0){
            break;
        }
        double cur_s = dp_speed_s[i];
        double cur_kappa = interpolate(std::vector ab(path_index2s.begin(), path_index2s.begin()+path_index2s_end_index), trajectory_init.points[1:path_index2s_end_index].kappa, cur_s)//插值找cur_s所对应的曲率
        auto max_speed = std::sqrt(max_lateral_accel / std::abs(cur_kappa) + 1e-10);
        double min_speed = 0;
        s_dot_lb[i] = min_speed;
        s_dot_ub[i] = max_speed;
    }

    for (int i=0; i<obs_st.s_in.size(); i++){
        auto obs_t = (obs_st.t_in[i] + obs_st.t_out[i]) / 2;    //取s t直线的中点， 作为obs_s obs_t的坐标
        auto obs_s = (obs_st.s_in[i] + obs_st.s_out[i]) / 2;
        auto obs_speed = (obs_st.s_out[i] - obs_st.s_in[i]) / (obs_st.t_out[i] - obs_st.t_in[i]);   //障碍物的纵向速度
        auto dp_s = interpolate([0, std::vector ab(dp_speed_t.begin(), dp_speed_t.begin()+dp_speed_end_index)], [0, std::vector ba(dp_speed_s.begin(), dp_speed_s.begin()+dp_speed_end_index), obs_t]);
        
        int j;
        for (j=0; j<dp_speed_t.size()-1; j++){
            if (dp_speed_t[0] > obs_st.t_in[i]) break;  //如果障碍物切入时间比0.5秒还要短，那么t_lb_index = 1
            else if (dp_speed_t[j] <= obs_st.t_in[i] && dp_speed_t[j+1])
            {   //否则遍历dp_speed_t 找到与obs_st_t_in_set[i]最近的点的编号
                break;
            }
        }
        int t_lb_index = j;

        for (j=0; j<dp_speed_t.size() - 1; j++){ //找到dp_speed_t中与obs_st_t_out_set[i]最近的时间，并将此时间的编号赋值给t_ub_index
            if dp_speed_t[0] > obs_st.t_out[i] break;
            else if (dp_speed_t[j] <= obs_st.t_out[i] && dp_speed_t[j+1] > obs_st.t_out[i])
            {
                break;
            }
        }
        int t_ub_index = j;
        t_lb_index = std::max(t_lb_index - 2, 3);   //最低为3， 因为碰瓷没法处理
        t_ub_index = std::min(t_ub_index + 2, dp_speed_end_index);

        if (obs_s > dp_s){
            for (int m=t_lb_index; m<=t_ub_index; m++){
                dp_t = dp_speed_t[m];
                s_ub[m] = std::min(s_ub[m], obs_st.s_in[i] + obs_speed * (dp_t - obs_st.t_in[i]);
            }
        }
        else{
            for (int m=t_lb_index; m<=t_ub_index; m++){
                dp_t = dp_speed_t[m];
                s_lb[m] = std::max(s_lb[m], obs_st.s_in[i] + obs_speed * (dp_t - obs_st.t_in[i]));
            }
        }
    }
}

void qp_velocity(){
    int dp_speed_end = 16;
    for (int i=0; i<dp_speed_s.size(); i++){
        if (dp_speed_s[i] < 0){
            dp_speed_end = i - 1;
        }
    }

    auto s_end = dp_speed_s[dp_speed_end];
    auto recommend_T = dp_speed_t[dp_speed_end];    //规划时间终点
    int qp_size = dp_speed_end + 1;  //qp的规模应该是dp的有效元素的个数 + 规划起点

    Eigen::MatrixXd Aeq = Eigen::Zeros(3*qp_size, 2*qp_size - 2);
    Eigen::MatrixXd beq = Eigen::Zeros(2*qp_size - 2, 1);
    Eigen::MatrixXd lb = Eigen::Ones(3*qp_size, 1);
    Eigen::MatrixXd ub = Eigen::Ones(3*qp_size, 1);
    Eigen::MatrixXd A = Eigen::Zeros(qp_size - 1, 3*qp_size);
    Eigen::MatrixXd b = Eigen::Zeros(qp_size - 1, 1);
    Eigen::MatrixXd H = Eigen::Zeros(3*n, 3*n);

    double dt = recommend_T / dp_speed_end;
    Eigen::MatrixXd A_sub(6, 2);
    double dt2 = std::pow(dt, 2);
    A_sub <<    1,         0,
                dt,         1,
                (1/3)*dt2,  (1/2)*dt,
                -1,         0,
                0,          -1,
                (1/6)*dt2,  dt/2;
    for (int i=0; i<qp_size-1; i++){
        Aeq(3*i:3*i+5, 2*i:2*i+2) = A_sub;
    }

    for (int i=0; i<qp_size-1; i++){
        A[i, 3*i] = 1;
        A[i, 3*i - 1] = -1;
    }

    for (int i=1; i<qp_size; i++){
        lb(3*i, 0) = s_lb[i];
        lb(3*i + 1, 0) = s_dot_lb[i];
        lb(3*i, 0) = -6;
        ub(3*i, 0) = s_ub[i];
        ub(3*i, 0) = s_dot_ub[i];
        ub(3*i, 0) = 4;
    }
    lb(0, 0) = 0;
    lb(1, 0) = start_point.frenet_s[1];
    lb(2, 0) = start_point.frenet_s[2];
    ub(0, 0) = lb(0, 0);
    ub(1, 0) = lb(1, 0);
    ub(2, 0) = lb(2, 0);

    Eigen::MatrixXd A_s_dot2 = Eigen::Zeros(3*qp_size, 3*qp_size);
    Eigen::MatrixXd A_jerk = Eigen::Zeros(3*qp_size, qp_size - 1);
    Eigen::MatrixXd A_ref = Eigen::Zeros(3*qp_size, 3*qp_size);
    Eigen::MatrixXd A4_sub(6, 1);
    A4_sub << 0, 0, 1, 0, 0, -1;
    
    for (int i=0; i<qp_size; i++){
        A_s_dot2(3*i+2, 3*i+2) = 1;
        A_ref(3*i+1, 3*i+1) = 1;
    }
    for (int i=0; i<qp_size-1; i++){
        A_jerk(3*i:3*i+5, i:i) = A4_sub;
    }

    auto H = w_cost_s_dot2 * A_s_dot2.dot(A_s_dot2.transpose()) + w_cost_jerk * A_jerk.dot(A_jerk.transpose()) + w_cost_v_ref * A_ref.dot(A_ref.transpose());
    H = 2 * H;

    Eigen::MatrixXd f = Eigen::Zeros(3*qp_size, 1);
    for (int i=0; i<qp_size; i++){
        f(3*i+1, 0) = -2 * w_cost_v_ref * speed_refernce
    }
    auto X = solve_qp(H, f, A, b, Aeq, beq, lb, ub, qp_size);
    // 缺乏曲率的约束，可以用|dl[i+1] - dl[i]|/ds <= kappa_max 近似约束曲率
    msgs::Trajectory trajectory_temp;
    for (int i=0; i<n; i++){
        trajectory_temp.points[i].frenet_s[0] = X[3*i];
        trajectory_temp.points[i].frenet_s[1] = X[3*i + 1];
        trajectory_temp.points[i].frenet_s[2] = X[3*i + 2];
        trajectory_temp.points[i].t = (i - 1) * dt;

    }
    add_density_s(trajectory_temp, trajectory_qp)
}

void save_msgs(msgs::TrajectoryPoint start, msgs::Object ego, msgs::ObjectList dynamic, msgs::Trajectory trajectory){
    start_point = start;
    ego_state = ego;
    dynamic_obstacles = dynamic;
    trajectory_init = trajectory;
}