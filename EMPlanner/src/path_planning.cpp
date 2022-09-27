#include "path_planning.h"
#include <map>
#include <cmath>
#include <eigen3/Eigen/Dense>
#include <msgs/ReferencePoint.h>
#include <msgs/TrajectoryPoint.h>
#include <msgs/Object.h>
#include <msgs/ObjectList.h>
#include <tools/tools.hpp>
#include <tools/frenet.hpp>
#include <string>
#include <ros/ros.h>
#include <vector>

int w_cost_collision, w_cost_smooth_dl, w_cost_smooth_ddl, w_cost_smooth_dddl;
int w_cost_ref, row, col;
double sample_s, smaple_l;
int w_cost_l, w_cost_dl, w_cost_ddl, w_cost_dddl, w_cost_centre;
int w_cost_end_l, w_cost_end_dl, w_cost_end_ddl;
int index2s[181] = {};

PathP::DP(){
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

void get_param(){
    
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

void DP::set_ego_state(const msgs::Object ego_data){
    ego_state = ego_data;
    // 处理自车信息
    msgs::ReferencePoint match_point, project_point;
    int match_point_index;
    auto item = obs_id2index.find(ego_state.id);
    if(item != obs_id2index.end()){
        int index = item->second;
        find_match_point(index-5, index+5, ego_state, rline, match_point, match_point_index, 5);
    }
    else{
        find_match_point(0, rline.points.size(), ego_state, rline, match_point, match_point_index, 5);
        obs_id2index.insert(std::make_pair(ego_state.id, match_point_index));
    }
    cal_project_point(ego_state, match_point, project_point);
    cartisian2frenet3D(ego_state, project_point);

    find_plan_start_point();
}

void DP::set_obstacles(const msgs::ObjectList obstacles_data){
    obstacles = obstacles_data;
    // 开始处理障碍物
    static_obstacles.objects.clear();
    dynamic_obstacles.objects.clear();
    static_obstacles.header = obstacles.header;
    dynamic_obstacles.header = obstacles.header;

    int n = 32;
    int count = 1;
    
    for (auto obs: obstacles.objects){
        double cos_theta = std::cos(ego_state.pose.theta), sin_theta = std::sin(ego_state.pose.theta);
        Eigen::Vector2d t_vector(cos_theta, sin_theta);     //自车heading的方向向量与法向量
        Eigen::Vector2d n_vector(-sin_theta, cos_theta);
        Eigen::Vector2d obs_vector(obs.pose.x - ego_state.pose.x, obs.pose.y - ego_state.pose.y);   //障碍物与自车的距离向量
        double lon_distance = obs_vector.dot(t_vector), lat_distance = obs_vector.dot(n_vector);    //障碍物与自车距离
       
        // 判断障碍是否在范围内，范围外障碍不考虑
        if(lon_distance < 60 && lon_distance > -10 && lat_distance < 10 && lat_distance > -10){
            // 判断障碍是否已经在上一周期查找过，如果查找过，则在上次匹配点前后继续查找，否则从头开始查找
            msgs::ReferencePoint match_point, project_point;
            int match_point_index;
            auto item = obs_id2index.find(obs.id);
            if(item != obs_id2index.end()){
                int index = item->second;
                find_match_point(index-10, index+10, obs, rline, match_point, match_point_index, 5);
            }
            else{
                find_match_point(0, rline.points.size(), obs, rline, match_point, match_point_index, 5);
                obs_id2index.insert(std::make_pair(obs.id, match_point_index));
            }
            // 根据匹配点计算投影点
            cal_project_point(obs, match_point, project_point);
            // 将障碍投影到frenet坐标下，坐标原点位参考线rline原点
            cartisian2frenet3D(obs, project_point);

            if(obs.frenet_s[1] < 0.1)
            static_obstacles.objects.push_back(obs);
            else
            dynamic_obstacles.objects.push_back(obs);
        }
    }
}

void calc_startpoint_stitch_trajectory(){
    ros::Time current_time = ros::Time::now();
    // 计算规划起点
    trajectory.line.clear();
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
        for(;index<pre_trajectory.line.size()-1;index++)
            if(pre_trajectory.line[index].header.stamp.toSec()<=current_time.toSec() && current_time.toSec()<pre_trajectory.line[index+1].header.stamp.toSec())
                break;
        
        double desired_x = pre_trajectory.line[index].pose.x;
        double desired_y = pre_trajectory.line[index].pose.y;
        double desired_theta = pre_trajectory.line[index].pose.theta;
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
            find_match_point(0, rline.points.size(), start_point, rline, match_point, match_point_index, 5);
            cal_project_point(start_point, match_point, project_point);
            cartisian2frenet3D(start_point, project_point);
            start_point.header.stamp = ros::Time().fromSec(current_time.toSec()+plan_dt);
        }
        else{      //控制正常
            int predict_index = index;
            for(;predict_index < pre_trajectory.line.size()-1; predict_index++)
                if(pre_trajectory.line[index].header.stamp.toSec()<=current_time.toSec()+plan_dt && current_time.toSec()+plan_dt<pre_trajectory.line[index+1].header.stamp.toSec())
                    break;
            start_point = pre_trajectory.line[predict_index];

            int length = 0;     //拼接轨迹
            for(int i=predict_index-1; i>=0; i--){
                trajectory.line.insert(trajectory.line.begin(), pre_trajectory.line[i]);
                length++;
                if(length > 20) break;
            }            
        }
    }

}



//     //计算index与s的转换，index2s（i）表示编号i对应的弧长s
// void index2s(msgs::ReferenceLine& path_pointline, msgs::ReferencePoint& origin_point, int& origin_match_point_index){
//     int n = 181;
//     int index2s[181] = {};

//     for (int i=1; i<=180; i++){
//         index2s[i] = std::sqrt(std::pow(path_pointline.points[i].x - path_pointline.points[i-1].x, 2) + std::pow(path_pointline[i].points.y - path_pointline[i-1].points.y, 2)) + index2s[i-1];
//     }
    
//     s0 = CalcSfromIndex2s(index2s, path_pointline, origin_pointline, origin_match_point_index);
//     for (int i=0; i<=180; i++){
//         index2s[i] = index2s[i] - s0;
//     }
//     return index2s;
    

// }
//     //计算点proj_x， proj_y对应的弧长，并判断投影点在匹配点的前面还是后面
// void CalcSfromIndex2s(int& index2s[181], msgs::ReferenceLine& path_pointline, msgs::ReferencePoint& project_point, int& proj_match_point_index){
//     Eigen::Vector2d vector1(project_point.x - path_pointline.points[proj_match_point_index].x, project_point.y - path_pointline.points[proj_match_point_index].y);
//     if (proj_match_point_index < path_pointline.points.size()){
//         Eigen::Vector2d vector2(path_pointline.points[proj_match_point_index + 1].x - path_pointline.points[proj_match_point_index].x, path_pointline.points[proj_match_point_index + 1].y - path_pointline.points[proj_match_point_index].y);
//     }
//     else{
//         Eigen::Vector2d vector2(path_pointline.points[proj_match_point_index].x - path_pointline.points[proj_match_point_index - 1].x, path_pointline.points[proj_match_point_index].y - path_pointline.points[proj_match_point_index - 1].y);
//     }
    
//     if (vector1.transpose().dot(vector2) > 0){      //投影点在匹配点前面
//         s = index2s[proj_match_point_index] + std::sqrt(vector1.transpose().dot(vector1));
//     }
//     else{
//         s = index2s[proj_match_point_index] - std::sqrt(vector1.transpose().dot(vector1));
//     }
// }



    //计算世界坐标系下的x_set， y_set上的点在frenet_path下的坐标s
void world2frenet_path(int& index2s[181], msgs::ReferenceLine& world_set, msgs::ReferenceLine& frenet_path_line, msgs::ReferenceLine& proj_set, msgs::ReferenceLine& project_point_set){
    int n = 128;
    int i;
    for (i=0; i<128; i++){
        if isnan(trajectory.points[i].frenet_l[0]){
            break;
        }
        trajectory.points[i].frenet_s[0] = CalcSFromIndex2S(index2s, frenet_path_line, proj_set.points[i], project_point_set.points[i]);
        Eigen::Vector2d n_r(-std::sin(proj_set.points[i].heading), std::cos(proj_set.points[i].heading));
        Eigen::Vector2d r_h(world_set.points[i].x - proj_set.points[i].x, world_set.points[i].y - proj_set.points[i].y);
        trajectory.points[i].frenet_l[0] = r_h.transpose().dot(n_r);
    }
}
    //计算frenet坐标系下的s_dot， l_dot(frenet_l[1]), dl/ds(frenet_l_dot[0])
void calc_dot_infrenet(msgs::Trajectory& trajectory, msgs::ReferenceLine& proj_set){
    int n = 128;
    for(int i=0; i<128; i++){
        if isnan(trajectory.points[i].frenet_l[0]){
            break;
        }
        Eigen::Vector2d v_h(trajectory.points[i].velocity.x, trajectory.points[i].velocity.y);
        Eigen::Vector2d n_r(-std::sin(proj_set.points[i].pose.theta), std::cos(proj_set.points[i].pose.theta));
        Eigen::Vector2d t_r(std::cos(proj_set.points[i].pose.theta), std::sin(proj_set.points[i].pose.theta));
        trajectroy.points[i].frenet_l[1] = v_h.transpose().dot(n_r);
        trajectroy.points[i].frenet_s[1] = v_h.transpose().dot(t_r/(1 - proj_set.points[i].kappa * trajectory.points[i].frenet_l[0]));
        if (std::abs(trajectory.points[i].frenet_s[1]) < 1e-6){
            trajectory.points[i].frenet_l_dot[0] = 0;
        }
        else{
            trajectory.points[i].frenet_l_dot[0] = trajectroy.points[i].frenet_l[1];
        }
    }
}
    //计算s_dot2, l_dot2, ddl
void calc_dot2_infrenet(msgs::Trajectory& trajectory, msgs::ReferenceLine& proj_set){
    int n = 128;
    for(int i=0; i<128; i++){
        if isnan(trajectory.points[i].frenet_l[0]){
            break;
        }
        Eigen::Vector2d a_h(trajectory.points[i].accel.x, trajectory.points[i].accel.y);
        Eigen::Vector2d n_r(-std::sin(proj_set.points[i].pose.theta), std::cos(proj_set.points[i].pose.theta));
        Eigen::Vector2d t_r(std::cos(proj_set.points[i].pose.theta), std::sin(proj_set.points[i].pose.theta));
        trajectroy.points[i].frenet_l[2] = a_h.transpose().dot(n_r - proj_set.points[i] * (1 - proj_set.points[i].kappa * trajectory.points[i].frenet_l[0]);
        trajectroy.points[i].frenet_s[2] = (1 / (1 - proj_set.points[i].pose.theta * trajectory.points[i].frenet_l[0])) * a_h.transpose().dot(t_r) + 2 * proj_set.points[i].kappa * trajectory.points[i].frenet_l_dot[0] * trajectory.points[i].frenet_s[2];
        if (std::abs(trajectory.points[i].frenet_s[1]) < 1e-6){
            trajectory.points[i].frenet_l_dot[1] = 0;
        }
        else{
            trajectory.points[i].frenet_l_dot[1] = (trajectroy.points[i].frenet_l[2] - trajectory.points[i].frenet_l_dot[0] * trajectory.points[i].frenet_s[2]) / std::pow(trajectory.points[i].frenet_s[1], 2);
        }
    }
}


void dp(){
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
        node_cost(i, 0) = CalcStartCost(start_point, cur_point, i, obs_set);
    }

    for(int j=1; j<col; j++){
        for (int i=0; i<row; i++){
            double cur_point.frenet_s[0] = start_point.frenet_s[0] + j * sample_s;
            double cur_point.frenet_l[0] = ((row + 1) / 2 - i) * smaple_l;
            for(int k=0; k<row; k++){
                double pre_point.frenet_s[0] = start_point.frenet_s + (j - 1) * sample_s;
                double pre_point.frenet_l[0] = ((row + 1) / 2 - k) * smaple_l;
                double cost_neighbour = CalcNeighbourCost(pre_point, cur_point, 999, obs_set);
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
    msg::Trajectory trajectory_temp;
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
    add_density(trajectory_temp, start_point)
}

double CalcCost(msgs::TrajectoryPoint& pre_point, msgs::TrajectoryPoint& cur_point, int cur_node, auto obs_set){
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
    auto l = a0 * cal + a1 * ds + a2 * std::pow(ds, 2) + a3 * std::pow(ds, 3) + a4 * std::pow(ds, 4) + a5 * std::pow(ds, 5);
    auto dl = a1 * cal + 2 * a2 * ds + 3 * a3 * std::pow(ds, 2) + 4 * a4 * std::pow(ds, 3) + 5 * a5 * std::pow(ds, 4);
    auto ddl = 2 * a2 * cal + 6 * a3 * ds + 12 * a4 * std::pow(ds, 2) + 20 * a5 * std::pow(ds, 3);
    auto dddl = 6 * a3 * cal + 25 * a4 * ds + 60 * a5 * std::pow(ds, 2);
 
    double cost_smooth = w_cost_smooth_dl * dl.transpose().dot(dl) + w_cost_smooth_ddl * ddl.transpose().dot(ddl) + w_cost_smooth_dddl * dddl.transpose().dot(dddl);
    double cost_ref = w_cost_ref * l.transpose().dot(l);
    double cost_collision = 0;
    
    for (int i=0; i<obs_set.rows(); i++){
        if (isnan(obs_set(i, 0))){
            break;
        }
        auto dlon = cal * obs_set(i, 0) - ds;
        auto dlat = cal * obs_set(i, 1) - l;
        auto square_d = std::pow(dlon, 2) + std::pow(dlat, 2);
        double cost_collision_once = CalcObsCost(square_d, w_cost_collision);
        cost_collision = cost_collision + cost_collision_once;
    }
    cost = cost_collision + cost_smooth + cost_ref;
    return cost;
}

void add_density(msgs::Trajectory& trajectory_temp, msgs::TrajectoryPoint& start_point){
    int ds = 1;
    std::vector s_cur, l_temp, dl_temp, ddl_temp;
    msgs::TrajectoryPoint end_point;
    msgs::TrajectoryPoint action_point;
    action_point = start_point
    for (int i=0; i<trajectory_temp.points.size(); i++){
        if (trajectory_temp.points[i].frenet_s[0] < 0){
            break;
        }
        for (int j=1; j<10000; j++){
            double s_node = action_point.frenet_s[0] + ds * (j - 1);
            if (s_node < trajectory_temp.points[i].frenet_s[0]){
                s_cur.push_back(s_node);
            }
            else{break;}
        }
        end_point.frenet_s[0] = trajectory_temp.points[i].frenet_s[0];
        end_point.frenet_l[0] = trajectory_temp.points[i].frenet_l[0];
        end_point.frenet_l[1] = 0;
        end_point.fernet_l[2] = 0;
        auto coeff = cal_quintic_coef(action_point, end_point);
        auto a0 = coeff[0], a1 = coeff[1], a2 = coeff[2], a3 = coeff[3], a4 = coeff[4], a5 = coeff[5];
        Eigen::Vector2d cal = Eigen::Ones(1, s_cur.size());
        auto l = a0 * cal + a1 * s_cur + a2 * std::pow(s_cur, 2) + a3 * std::pow(s_cur, 3) + a4 * std::pow(s_cur, 4) + a5 * std::pow(s_cur, 5);
        auto dl = a1 * cal + 2 * a2 * s_cur + 3 * a3 * std::pow(s_cur, 2) + 4 * a4 * std::pow(s_cur, 3) + 5 * a5 * std::pow(s_cur, 4);
        auto ddl = 2 * a2 * cal + 6 * a3 * s_cur + 12 * a4 * std::pow(s_cur, 2) + 20 * a5 * std::pow(s_cur, 3);
        l_temp.push_back(l);
        dl_temp.push_back(dl);
        ddl_temp.push_back(ddl);
        action_point = end_point;
        s_cur.clear();
    }
    for (int i=0; i<s_cur.size(); i++){
        if (i == 60){
            break;
        }
        trajectory.points[i].frenet_s[0] = s_cur[i];
        trajectory.points[i].frenet_l[0] = l_temp[i];
        trajectory.points[i].frenet_l[1] = dl_temp[i];
        trajectory.points[i].frenet_l[2] = ddl_temp[i];
    }
}

void path_frenet2cartesian(msg::Trajectory& trajectory, msgs::ReferenceLine& rline){
    msgs::Object cal_object;
    msgs::ReferencePoint proj_point;
    msgs::ReferencePoint match_point;
    for (int i=0; i<trajectory.points.size(); i++){
        if (isnan(trajectory.point[i].frenet_s[0])){
            break;
        }
        cal_object.location = trajectory.points[i];
        find_match_point(0, rline.points.size(), cal_object, rline, match_point, match_point_index, 10);
        cal_project_point_f2c(cal_object, match_point, proj_point);
        trajectory_ad.points[i].pose.x = proj_point.x + trajectory.points[i].frenet_l[0] * -std::sin(proj_point.theta);
        trajectory_ad.points[i].pose.y = proj_point.y + trajectory.points[i].frenet_l[0] * std::cos(proj_point.theta);
        trajectory_ad.points[i].pose.theta = proj_point.theta + std::atan(trajectory.points[i].frenet_l[1] / (1 - proj_point.kappa * trajectory.points[i].frenet_l[0]));
        trajectory.points[i].kappa = ((trajectory.points[i].frenet_l[2] + proj_point.theta * trajectory.points[i].frenet_l[1] * std::tan(trajectory.points[i].pose.theta - proj_point.theta)) * std::pow(std::cos(trajectory.points[i].pose.theta - proj_point.theta), 2) / (1 - proj_point.kappa * trajectory.points[i].frenet_l[0]) + proj_point.kappa) * std::cos(trajectory.points[i].pose.theta - proj_point.theta) / (1 - proj_point.kappa * trajectory.points[i].frenet_l[0]);
    }
}

////
//qp
////
Eigen::Vector2d l_min = Eigen::Ones(n ,1);
Eigen::Vector2d l_max = Eigen::Ones(n ,1);

void get_boundary(){
    int n = 60;
    l_min = l_min * -6;
    l_max = l_max * 6;

    for (auto obs: static_obstacles.objects){
        double obs_s_min = obs.location.frenet_s[0] - obs.length / 2;
        double obs_s_max = obs.location.frenet_s[0] + obs.length / 2;
        int start_index = find_near_index(trajectory, obs_s_min);
        int end_index = find_near_index(trajectory, obs_s_max);
        int centre_index = find_near_index(trajectory, obs.location.frenet_s[0]);

        if (start_index == 1 && end_index == 1){
            continue;
        }
        double path_l = trajectory.points[centre_index].frenet_l[0];
        if (path_l > obs.location.frenet_l[0]){
            for (int j=start_index; j<=end_index; j++){
                l_min(j, 0) = std::max(l_min(j, 0), l_min(j, 0) + obs.location.frenet_l[0] + obs.width/2);
            }
        }
        else{
            for (int j=start_index; j<=end_index; j++){
                l_max(j, 0) = std::min(l_max(j, 0), obs.location.frenet_l[0] - obs.width/2);

            }
        }
    }
}

int find_near_index(msg::Trajectory& trajectory, double obs_s){
    int y = 0;
    int index = 0;
    if (trajectroy.points[0].frenet_s[0] > = obs_s){
        y = 0;
    }
    else if (trajectory.points[-1].frenet_s[0] < obs_s)
    {
        y = 60;
    }
    else{
        while (trajectory.points[index] < obs_s){
            index = index + 1;
        }
    }

    if (trajectory.points[index] - obs_s > obs_s - trajectory.points[index - 1]){
        y = index - 1;
    }
    else{ y = index; }
    return y;
    
}

void qp(){
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
        row_temp = 2*i - 1 + 2;
        col_temp = 3*i - 2 + 3;
        Aeq(row_temp: row_temp + 2 , col_temp: col_temp+6) = Aeq_sub;
    }
    for (int i=0, i<n; i++){
    // for (int i=1, i<n; i++){     //二次规划解决办法之二
        row_temp = 8*i - 7 + 8;
        col_temp = 3*i - 2 + 3;
        A(row_temp: row_temp + 8 , col_temp: col_temp+3) = A_sub;
    }
    int front_index = std::ceil(d1/ds);
    int back_index = std::ceil(d2/ds);

    for (int i=0; i<n; i++){
    // for (int i=1, i<n; i++){     //二次规划解决办法之二
        int index1 = std::min(i + front_index, n);
        int index2 = std::max(i - back_index, 1);
        b(8*i - 7:  8*i + 1, 0) << l_max(index1, 0) - w/2, l_max(index1, 0) + w/2, l_max(index2, 0) - w/2, l_max(index2, 0) + w/2, -l_min(index1, 0) + w/2, -l_min(index1, 0) - w/2, -l_min(index2, 0) + w/2, -l_min(index2, 0) - w/2;
    }
    Eigen::MatrixXd lb = Eigen::Ones(3*n, 1);
    Eigen::MatrixXd ub = Eigen::Ones(3*n, 1);
    lb *= -99999;
    ub *= 99999;
    lb(0, 0) = start_point.frenet_l[9];
    lb(1, 0) = start_point.frenet_l[1];
    lb(2, 0) = start_point.frenet_l[2];
    ub = lb;

    for (int i=0; i<n; i++){
        H_L(3*i - 2, 3*i - 2) = 1;
        H_DL(3*i - 1, 3*i - 1) = 1;
        H_DDL(3*i, 3*i) = 1;
    }

    Eigen::MatrixXd H_CENTRE = H_L;
    Eigen::VectorXd H_dddl_sub(6) << 0, 0, 1, 0, 0, -1;
    for (int i=0; i<n-1; i++){
        row_temp = i;
        col_temp = 3*i - 2;
        H_DDDL(row_temp, col_temp: col_temp + 6) = H_dddl_sub;
    }

    H_L_END(3*n - 2, 3*n - 2) = 1;
    H_DL_END(3*n - 1, 3*n - 1) = 1;
    H_DDL_END(3*n, 3*n) = 1;
    H = w_cost_l * H_L.transpose().dot(H_L) + w_cost_dl * H_DL.transpose().dot(H_DL) + w_cots_ddl * H_DDL.transpose().dot(H_DDL) + w_cost_dddl * H_DDDL.transpose().dot(H_DDDL) + w_cost_centre * H_CENTRE.transpose().dot(H_CENTRE) + w_cost_end_l * H_L_END.transpose().dot(H_L_END) + w_cost_end_dl * H_DL_END.transpose().dot(H_DL_END) + w_cost_end_ddl * H_DDL_END.transpose().dot(H_DDL_END);
    H = 2 * H;

    Eigen::VectorXd f = Eigen::Zeros(3*n);
    auto centre_line = 0.5 * (l_min + l_max);
    // auto centre_line = trajectory.frenet_l[0];       //二次规划奔溃解决方法之一
    
    for (int i=0; i<n; i++){
        f[3*i - 2] = -2 * centre_line[i];
    }
    f = w_cost_centre * f;

    f[3*n - 2] -= 2 * end_l_desire * w_cost_end_l;
    f[3*n - 1] -= 2 * end_dl_desire * w_cost_end_dl;
    f[3*n] -= 2 * end_ddl_desire * w_cost_end_ddl;

    X = solve(H, f, A, b, Aeq, beq, lb, ub);
    // 缺乏曲率的约束，可以用|dl[i+1] - dl[i]|/ds <= kappa_max 近似约束曲率
    for (int i=0; i<0; i++){
        trajectory_qp.points[i].frenet_l[0] = X[3*i];
        trajectory_qp.points[i].frenet_l[1] = X[3*i + 1];
        trajectory_qp.points[i].frenet_l[2] = X[3*i + 2];
    }
}
