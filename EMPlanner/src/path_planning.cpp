#include "path_planning.h"
#include <map>
#include <cmath>
#include <eigen3/Eigen/Dense>
#include <string>
#include <ros/ros.h>
#include <vector>
#include "OsqpEigen/OsqpEigen.h"
#include <msgs/ReferencePoint.h>
#include <msgs/ReferenceLine.h>
#include <msgs/TrajectoryPoint.h>
#include <msgs/Object.h>
#include <msgs/ObjectList.h>
#include <tools/tools.hpp>
#include <tools/frenet.hpp>
#include <tools/pathp_tools.hpp>

int w_cost_collision, w_cost_smooth_dl, w_cost_smooth_ddl, w_cost_smooth_dddl;
int w_cost_ref, row, col;
double sample_s, smaple_l;
int w_cost_l, w_cost_dl, w_cost_ddl, w_cost_dddl, w_cost_centre;
int w_cost_end_l, w_cost_end_dl, w_cost_end_ddl;

PathP::PathP(){
    first_run = true;
}

void PathP::get_referline(const msgs::ReferenceLine rline_data){
    rline = rline_data;
}

void PathP::set_dt(double planning_delta_time){
    plan_dt = planning_delta_time;
}

void PathP::get_trajactory(const msgs::Trajectory trajectory_data){
    trajectory = trajectory_data;
}

void PathP::get_param(){
    // if (!ros::param::get("w_cost_collision", w_cost_collision))
    //     ROS_WARN_STREAM("Cannot load w_cost_collision. Standard value is: " << w_cost_collision);
    ros::param::get("/row", row);
    ros::param::get("/col", col);
    ros::param::get("/sample_s", sample_s);
    ros::param::get("/smaple_l", smaple_l);
    ros::param::get("/w_cost_collision", w_cost_collision);
    ros::param::get("/w_cost_smooth_dl", w_cost_smooth_dl);
    ros::param::get("/w_cost_smooth_ddl", w_cost_smooth_ddl);
    ros::param::get("/w_cost_smooth_dddl", w_cost_smooth_dddl);
    ros::param::get("/w_cost_ref", w_cost_ref);
    ros::param::get("/w_cost_l", w_cost_l);
    ros::param::get("/w_cost_dl", w_cost_dl);
    ros::param::get("/w_cost_ddl", w_cost_ddl);
    ros::param::get("/w_cost_dddl", w_cost_dddl);
    ros::param::get("/w_cost_centre", w_cost_centre);
    ros::param::get("/w_cost_end_l", w_cost_end_l);
    ros::param::get("/w_cost_end_dl", w_cost_end_dl);
    ros::param::get("/w_cost_end_ddl", w_cost_end_ddl);
}

void PathP::path_planning(const msgs::Object ego_data, const msgs::ObjectList obstacles_data){
    get_param();
    set_ego_state(*ego_data);
    calc_startpoint_stitch_trajectory();
    set_obstacles(*obstacles_data);
    dp();
    get_boundary();
    qp();
    trajectory = trajectory_qp;
    path_frenet2cartesian(trajectory, rline, index2s);
    pre_trajectory = trajectory;
}

void PathP::set_ego_state(const msgs::Object ego_data){
    ego_state = ego_data;
    // 处理自车信息
    msgs::ReferencePoint match_point, project_point;
    int match_point_index;
    auto item = obs_id2index.find(ego_state.id);
    if(item != obs_id2index.end()){
        int index = item->second;
        cal_match_point(index-5, index+5, ego_state, rline, match_point, match_point_index, 5);
    }
    else{
        cal_match_point(0, rline.points.size(), ego_state, rline, match_point, match_point_index, 5);
        obs_id2index.insert(std::make_pair(ego_state.id, match_point_index));
    }
    cal_project_point(ego_state, match_point, project_point);
    index2s(match_point_index);
    world2frenet_path(index2s, ego_state, rilne, project_point, match_point_index, 2);
}

void PathP::set_obstacles(const msgs::ObjectList obstacles_data){
    obstacles = obstacles_data;
    // 开始处理障碍物
    static_obstacles.objects.clear();
    dynamic_obstacles.objects.clear();
    static_obstacles.header = obstacles.header;
    dynamic_obstacles.header = obstacles.header;

    for (auto obs: obstacles.objects){
        double cos_theta = std::cos(ego_state.pose.theta), sin_theta = std::sin(ego_state.pose.theta);
        Eigen::Vector2d t_vector(cos_theta, sin_theta);     //自车heading的方向向量与法向量
        Eigen::Vector2d n_vector(-sin_theta, cos_theta);
        Eigen::Vector2d obs_vector(obs.pose.x - ego_state.pose.x, obs.pose.y - ego_state.pose.y);   //障碍物与自车的距离向量
        double lon_distance = obs_vector.dot(t_vector), lat_distance = obs_vector.dot(n_vector);    //障碍物与自车距离
       
        // 判断障碍是否在范围内，范围外障碍不考虑
        if(lon_distance < 60 && lon_distance > -10 && lat_distance < 25 && lat_distance > -25){ //之前的lat值为±10，对于动态障碍物来说太小了
            // 判断障碍是否已经在上一周期查找过，如果查找过，则在上次匹配点前后继续查找，否则从头开始查找
            msgs::ReferencePoint match_point, project_point;
            int match_point_index;
            auto item = obs_id2index.find(obs.id);
            if(item != obs_id2index.end()){
                int index = item->second;
                cal_match_point(index-10, index+10, obs, rline, match_point, match_point_index, 5);
            }
            else{
                cal_match_point(0, rline.points.size(), obs, rline, match_point, match_point_index, 5);
                obs_id2index.insert(std::make_pair(obs.id, match_point_index));
            }
            // 根据匹配点计算投影点
            cal_project_point(obs, match_point, project_point);
            // 将障碍投影到frenet坐标下，坐标原点位参考线rline原点
            world2frenet_path(index2s, obs, rilne, project_point, match_point_index, 2);

            if(obs.frenet_s[1] < 0.1)
            static_obstacles.objects.push_back(obs);
            else
            dynamic_obstacles.objects.push_back(obs);
        }
    }
}

void PathP::calc_startpoint_stitch_trajectory(){
    ros::Time current_time = ros::Time::now();
    // 计算规划起点
    trajectory.points.clear();
    if(first_run){
        start_point.pose = ego_state.pose;
        start_point.frenet_s = ego_state.frenet_s;
        start_point.frenet_d = ego_state.frenet_d;
        start_point.frenet_d_dot = ego_state.frenet_d_dot;
        start_point.header.stamp = ros::Time().fromSec(current_time.toSec()+plan_dt);
        first_run = false;
    }
    else{
        int index=0;    //上一周期的车辆理论位置
        for(;index<pre_trajectory.points.size()-1;index++)
            if(pre_trajectory.points[index].header.stamp.toSec()<=current_time.toSec() && current_time.toSec()<pre_trajectory.points[index+1].header.stamp.toSec())
                break;
        
        double desired_x = pre_trajectory.points[index].pose.x;
        double desired_y = pre_trajectory.points[index].pose.y;
        double desired_theta = pre_trajectory.points[index].pose.theta;
        Eigen::Vector2d t_vector(std::cos(desired_theta), std::sin(desired_theta));
        Eigen::Vector2d n_vector(-std::sin(desired_theta), std::cos(desired_theta));
        Eigen::Vector2d err_vector(ego_state.pose.x-desired_x, ego_state.pose.y-desired_y);
        double lon_err = std::abs(err_vector.transpose().dot(t_vector));
        double lat_err = std::abs(err_vector.transpose().dot(n_vector));

        // //当车辆速度加速度不是世界坐标系时
        // ego_state.velocity.x = ego_state.velocity.x * std::cos(ego_state.pose.theta) - ego_state.velocity.y * std::sin(ego_state.pose.theta)
        // ego_state.velocity.y = ego_state.velocity.x * std::sin(ego_state.pose.theta) + ego_state.velocity.y * std::cos(ego_state.pose.theta)
        // ego_state.accel.x = ego_state.accel.x * std::cos(ego_state.pose.theta) - ego_state.accel.y * std::sin(ego_state.pose.theta)
        // ego_state.accel.y = ego_state.accel.x * std::cos(ego_state.pose.theta) + ego_state.accel.y * std::sin(ego_state.pose.theta)

        if(lon_err > 2.5 || lat_err > 0.5){     //误差过大，动力学外推
            start_point.pose.x = ego_state.pose.x + ego_state.velocity.x * plan_dt + 0.5 * ego_state.accel.x * plan_dt * plan_dt;
            start_point.pose.y = ego_state.pose.y + ego_state.velocity.y * plan_dt + 0.5 * ego_state.accel.y * plan_dt * plan_dt;
            start_point.pose.theta = ego_state.pose.theta + ego_state.velocity.theta * plan_dt + 0.5 * ego_state.accel.theta * plan_dt * plan_dt;
            start_point.velocity.x = ego_state.velocity.x + ego_state.accel.x * plan_dt;
            start_point.velocity.y = ego_state.velocity.y + ego_state.accel.y * plan_dt;
            start_point.velocity.theta = ego_state.velocity.theta + ego_state.accel.theta * plan_dt;
            start_point.accel = ego_state.accel;
            start_point.kappa = 0;

            msgs::ReferencePoint match_point, project_point;
            int match_point_index;
            cal_match_point(0, rline.points.size(), start_point, rline, match_point, match_point_index, 5);
            cal_project_point(start_point, match_point, project_point);
            world2frenet_path(index2s, start_point, rilne, project_point, match_point_index, 3);
            start_point.header.stamp = ros::Time().fromSec(current_time.toSec()+plan_dt);
        }
        else{      //控制正常
            int predict_index = index;
            for(;predict_index < pre_trajectory.points.size()-1; predict_index++)
                if(pre_trajectory.points[index].header.stamp.toSec()<=current_time.toSec()+plan_dt && current_time.toSec()+plan_dt<pre_trajectory.points[index+1].header.stamp.toSec())
                    break;
            start_point = pre_trajectory.points[predict_index];

            int length = 0;     //拼接轨迹
            for(int i=predict_index-1; i>=0; i--){
                trajectory.points.insert(trajectory.points.begin(), pre_trajectory.points[i]);
                length++;
                if(length > 20) break;
            }            
        }
    }
}

    //计算index与s的转换，index2s（i）表示编号i对应的弧长s
void PathP::index2s(int& origin_match_point_index){
    int n = 180;
    index2s[0] = 0;
    msgs::ReferencePoint origin_point;
    origin_point.x = ego_state.pose.x;
    origin_point.y = ego_state.pose.y;

    for (int i=1; i<=n; i++){
        index2s[i] = std::sqrt(std::pow(rilne.points[i].x - rilne.points[i-1].x, 2) + std::pow(rilne[i].points.y - rilne[i-1].points.y, 2)) + index2s[i-1];
    }
    s0 = CalcSfromIndex2s(index2s, rilne, origin_point, origin_match_point_index);
    for (int i=0; i<=n; i++){
        index2s[i] = index2s[i] - s0;
    }
}

////
//dp
////

void PathP::dp(){   //输出的s， l不包含规划起点
    Eigen::MatrixXd <double, row, col> node_cost;
    Eigen::MatrixXd <double, row, col> pre_node_index;
    node_cost.setOnes();
    node_cost = node_cost * 99999;

    Eigen::Vector2d obs_set = Eigen::Zero(static_obstacles.objects.size()+dynamic_obstacles.objects.size(), 2);
    for (int i=0; i<static_obstacles.objects.size(); i++){
        obs_set(i, 0) = static_obstacles.objects[i].frenet_s[0];
        obs_set(i, 1) = static_obstacles.objects[i].frenet_l[0];
    }
    for (int i=0; i<dynamic_obstacles.objects.size(); i++){
        obs_set(i+static_obstacles.objects.size(), 0) = dynamic_obstacles.objects[i].frenet_s[0];
        obs_set(i+static_obstacles.objects.size(), 1) = dynamic_obstacles.objects[i].frenet_l[0];
    }

    msgs::TrajectoryPoint pre_point;
    msgs::TrajectoryPoint cur_point;
    for (int i=0; i<row; i++){
        node_cost(i, 0) = PathP::CalcCost(start_point, cur_point, i, obs_set);
    }

    for(int j=1; j<col; j++){
        for (int i=0; i<row; i++){
            double cur_point.frenet_s[0] = start_point.frenet_s[0] + j * sample_s;
            double cur_point.frenet_l[0] = ((row + 1) / 2 - i) * smaple_l;
            for(int k=0; k<row; k++){
                double pre_point.frenet_s[0] = start_point.frenet_s + (j - 1) * sample_s;
                double pre_point.frenet_l[0] = ((row + 1) / 2 - k) * smaple_l;
                double cost_neighbour = PathP::CalcCost(pre_point, cur_point, 999, obs_set);
                double pre_min_cost = node_cost(k, j-1);
                if ((pre_min_cost + cost_neighbour) < node_cost(i, j)){
                    node_cost(i, j) = pre_min_cost + cost_neighbour;
                    pre_node_index(i, j) = k;
                }
            } 
        }
    }
    int index = 0;
    int min_cost = 99999;
    msgs::Trajectory trajectory_temp;
    for(int i=0; i<row; i++){
        if (node_cost(i, col-1) < min_cost){
            min_cost = node_cost(i, col-1);
            index = i;
        }

        int dp_node_list_row[col] = {};
        for(int i=0; i<col; i++){
            int pre_index = pre_node_index(index, col-i);
            dp_node_list_row[col - i] = index;
            index = pre_index;
        }

        for (int i=0; i<col; i++){
            trajectory_temp.points[i].frenet_s[0] = start_point.frenet_s[0] + i * sample_s;
            trajectory_temp.points[i].frenet_l[0] = ((row + 1) / 2 - dp_node_list_row[i]) * smaple_l;
        }
    }
    add_density_dp(trajectory_temp, start_point, trajectory_dp)
}

double CalcCost(auto& pre_point, auto& cur_point, int cur_node, auto obs_set){
    Eigen::Vector2d ds = Eigen::Zero(row+1, 1), l = Eigen::Zero(row+1, 1), dl = Eigen::Zero(row+1, 1);
    Eigen::Vector2d ddl = Eigen::Zero(row+1, 1), dddl = Eigen::Zero(row+1, 1), cal = Eigen::Ones(row+1, 1);
    
    if (cur_node == 999){
        pre_point.frenet_l[1] = 0;
        pre_point.frenet_l[2] = 0;
        cur_point.frenet_l[1] = 0;
        cur_point.frenet_l[2] = 0;
        for (int i=0; i<row+1; i++){
            ds(i, 0) = pre_point.frenet_s[0] + (i - 1) * (cur_point.frenet_s[0] - pre_point.frenet_s[0]) / 10;
        }
    }
    else{
        cur_point.frenet_l[0] = ((row + 1) / 2 - cur_node) * sample_l;
        cur_point.frenet_l[1] = 0;
        cur_point.frenet_l[2] = 0;
        cur_point.frenet_s[0] = pre_point.frenet_s[0] + sample_s;
        for (int i=0; i<row+1; i++){
            ds(i, 0) = pre_point.frenet_s[0] + (i - 1) * sample_s / 10;
        }
    }
    auto coeff = cal_quintic_coef(pre_point, cur_point);
    auto a0 = coeff[0], a1 = coeff[1], a2 = coeff[2], a3 = coeff[3], a4 = coeff[4], a5 = coeff[5];
    auto l = a0 * cal + a1 * ds + a2 * cal_pow(ds, 2) + a3 * cal_pow(ds, 3) + a4 * cal_pow(ds, 4) + a5 * cal_pow(ds, 5);
    auto dl = a1 * cal + 2 * a2 * ds + 3 * a3 * cal_pow(ds, 2) + 4 * a4 * cal_pow(ds, 3) + 5 * a5 * cal_pow(ds, 4);
    auto ddl = 2 * a2 * cal + 6 * a3 * ds + 12 * a4 * cal_pow(ds, 2) + 20 * a5 * cal_pow(ds, 3);
    auto dddl = 6 * a3 * cal + 25 * a4 * ds + 60 * a5 * cal_pow(ds, 2);
 
    double cost_smooth = w_cost_smooth_dl * dl.transpose().dot(dl) + w_cost_smooth_ddl * ddl.transpose().dot(ddl) + w_cost_smooth_dddl * dddl.transpose().dot(dddl);
    double cost_ref = w_cost_ref * l.transpose().dot(l);
    double cost_collision = 0;
    
    for (int i=0; i<obs_set.rows(); i++){
        auto dlon = cal * obs_set(i, 0) - ds;
        auto dlat = cal * obs_set(i, 1) - l;
        auto square_d = cal_pow(dlon, 2) + cal_pow(dlat, 2);
        double cost_collision_once = CalcObsCost(square_d, w_cost_collision);
        cost_collision = cost_collision + cost_collision_once;
    }
    cost = cost_collision + cost_smooth + cost_ref;
    return cost;
}


////
//qp
////
std::vector l_min(60, -6);
std::vector l_max(60, 6);

void PathP::get_boundary(){
    int n = 60;
    for (int i=0; i<n; i++){
        l_min[i] = -6;
        l_max[i] = 6;
    }

    for (auto obs: static_obstacles.objects){
        double obs_s_min = obs.frenet_s[0] - obs.length / 2;
        double obs_s_max = obs.frenet_s[0] + obs.length / 2;
        int start_index = find_near_index(trajectory_dp, obs_s_min);
        int end_index = find_near_index(trajectory_dp, obs_s_max);
        int centre_index = find_near_index(trajectory_dp, obs.frenet_s[0]);

        if (start_index == 1 && end_index == 1){
            continue;
        }
        double path_l = trajectory_dp.points[centre_index].frenet_l[0];
        if (path_l > obs.frenet_l[0]){
            for (int j=start_index; j<=end_index; j++){
                l_min[j] = std::max(l_min[j], l_min[j] + obs.frenet_l[0] + obs.width/2);
            }
        }
        else{
            for (int j=start_index; j<=end_index; j++){
                l_max[j] = std::min(l_max[j], obs.frenet_l[0] - obs.width/2);
            }
        }
    }


void PathP::qp(){
    int n = 60;
    Eigen::MatrixXd H_L = Eigen::Zeros(3*n, 3*n);
    Eigen::MatrixXd H_DL = Eigen::Zeros(3*n, 3*n);
    Eigen::MatrixXd H_DDL = Eigen::Zeros(3*n, 3*n);
    Eigen::MatrixXd H_DDDL = Eigen::Zeros(n-1, 3*n);
    Eigen::MatrixXd H_L_END = Eigen::Zeros(3*n, 3*n);
    Eigen::MatrixXd H_DL_END = Eigen::Zeros(3*n, 3*n);
    Eigen::MatrixXd H_DDL_END = Eigen::Zeros(3*n, 3*n);
    Eigen::MatrixXd Aeq = Eigen::Zeros(2*n-2, 3*n);
    Eigen::MatrixXd beq = Eigen::Zeros(2*n-2, 1);
    Eigen::MatrixXd A = Eigen::Zeros(8*n, 3*n);
    Eigen::MatrixXd b = Eigen::Zeros(8*n, 1);
    Eigen::MatrixXd H = Eigen::Zeros(3*n, 3*n)
    int end_l_desire = 0, end_dl_desire = 0, end_ddl_desire = 0;

    int ds = 1;
    Eigen::MatrixXd Aeq_sub(2, 6);
    Eigen::MatrixXd A_sub(8, 3);
    Aeq_sub <<  1, ds,  std::pow(ds, 2/3),  -1, 0,  std::pow(ds, 1/3), 
                0, 1,   ds/2,               0,  -1, ds/2;     
    
    double d1 = ego_state.length / 2, d2 = ego_state.length / 2, w = ego_state.width;
    
    A_sub <<    1,  d1,     0,
                1,  d1,     0,
                1,  -d2,    0,
                1,  -d2,    0,
                -1, -d1,    0,
                -1, -d1,    0,
                -1, d2,     0, 
                -1, d2,     0;
    int row_temp, col_temp;
    for (int i=0, i<n-1; i++){
        row_temp = 2*i;
        col_temp = 3*i;
        Aeq(row_temp: row_temp + 2 , col_temp: col_temp+6) = Aeq_sub;
    }
    for (int i=0, i<n; i++){
    // for (int i=1, i<n; i++){     //二次规划解决办法之二
        row_temp = 8*i;
        col_temp = 3*i;
        A(row_temp: row_temp + 8 , col_temp: col_temp+3) = A_sub;
    } 
    int front_index = std::ceil(d1/ds);
    int back_index = std::ceil(d2/ds);

    for (int i=0; i<n; i++){
    // for (int i=1, i<n; i++){     //二次规划解决办法之二
        int index1 = std::min(i + front_index, n);
        int index2 = std::max(i - back_index, 1);
        b(8*i:  8*i + 8, 0) << l_max[index1] - w/2, l_max[index1] + w/2, l_max[index2] - w/2, l_max[index2] + w/2, -l_min[index1] + w/2, -l_min[index1] - w/2, -l_min[index2] + w/2, -l_min[index2] - w/2;
    }
    Eigen::MatrixXd lb = Eigen::Ones(3*n, 1);
    Eigen::MatrixXd ub = Eigen::Ones(3*n, 1);
    lb *= -99999;
    ub *= 99999;
    lb(0, 0) = start_point.frenet_l[0];
    lb(1, 0) = start_point.frenet_l[1];
    lb(2, 0) = start_point.frenet_l[2];
    ub = lb;

    for (int i=0; i<n; i++){
        H_L(3*i, 3*i) = 1;
        H_DL(3*i + 1, 3*i + 1) = 1;
        H_DDL(3*i + 2, 3*i + 2) = 1;
    }

    Eigen::MatrixXd H_CENTRE = H_L;
    Eigen::VectorXd H_dddl_sub(6) << 0, 0, 1, 0, 0, -1;
    for (int i=0; i<n-1; i++){
        row_temp = i;
        col_temp = 3*i;
        H_DDDL(row_temp, col_temp: col_temp + 6) = H_dddl_sub;
    }

    H_L_END(3*n - 2, 3*n - 2) = 1;
    H_DL_END(3*n - 1, 3*n - 1) = 1;
    H_DDL_END(3*n, 3*n) = 1;
    H = w_cost_l * H_L.transpose().dot(H_L) + w_cost_dl * H_DL.transpose().dot(H_DL) + w_cots_ddl * H_DDL.transpose().dot(H_DDL) + w_cost_dddl * H_DDDL.transpose().dot(H_DDDL) \
        + w_cost_centre * H_CENTRE.transpose().dot(H_CENTRE) + w_cost_end_l * H_L_END.transpose().dot(H_L_END) + w_cost_end_dl * H_DL_END.transpose().dot(H_DL_END) + w_cost_end_ddl * H_DDL_END.transpose().dot(H_DDL_END);
    H = 2 * H;

    Eigen::VectorXd f = Eigen::Zeros(3*n);
    std::vector centre_line;
    for (int i=0; i<l_min.size(); i++){
        centre_line[i] = 0.5 * (l_min[i] + l_max[i]);
    }
    // auto centre_line = trajectory_dp.frenet_l[0];       //二次规划奔溃解决方法之一
    
    for (int i=0; i<n; i++){
        f[3*i] = -2 * centre_line[i]; 
    }
    f = w_cost_centre * f;

    f[3*n - 2] -= 2 * end_l_desire * w_cost_end_l;
    f[3*n - 1] -= 2 * end_dl_desire * w_cost_end_dl;
    f[3*n] -= 2 * end_ddl_desire * w_cost_end_ddl;

    auto X = solve_qp(H, f, A, b, Aeq, beq, lb, ub, n);
    // 缺乏曲率的约束，可以用|dl[i+1] - dl[i]|/ds <= kappa_max 近似约束曲率
    msgs::Trajectory trajectory_temp;
    for (int i=0; i<n; i++){
        trajectory_temp.points[i].frenet_l[0] = X[3*i];
        trajectory_temp.points[i].frenet_l[1] = X[3*i + 1];
        trajectory_temp.points[i].frenet_l[2] = X[3*i + 2];
    }
    add_density_qp(trajectory_temp, trajectory_qp);
}

void get_msgs(){
    return start_point, ego_state, dynamic_obstacles, trajectory
}
