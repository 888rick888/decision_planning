#pragma once
#include <cmath>
#include <eigen3/Eigen/Dense>
#include "OsqpEigen/OsqpEigen.h"
#include <msgs/Object.h>
#include <msgs/ReferenceLine.h>
#include <msgs/ReferencePoint.h>
#include <tools.hpp>
#include <vector>


    //计算点proj_x， proj_y对应的弧长，并判断投影点在匹配点的前面还是后面
template <class T>
void CalcSfromIndex2s(std::vector<double>& index2s, msgs::ReferenceLine& rilne, msgs::ReferencePoint& project_point, int& proj_match_point_index){
    Eigen::Vector2d vector1(project_point.x - rilne.points[proj_match_point_index].x, project_point.y - rilne.points[proj_match_point_index].y);
    if (proj_match_point_index < rilne.points.size()){
        Eigen::Vector2d vector2(rilne.points[proj_match_point_index + 1].x - rilne.points[proj_match_point_index].x, rilne.points[proj_match_point_index + 1].y - rilne.points[proj_match_point_index].y);
    }
    else{
        Eigen::Vector2d vector2(rilne.points[proj_match_point_index].x - rilne.points[proj_match_point_index - 1].x, rilne.points[proj_match_point_index].y - rilne.points[proj_match_point_index - 1].y);
    }
    if (vector1.transpose().dot(vector2) > 0){      //投影点在匹配点前面
        s = index2s[proj_match_point_index] + std::sqrt(vector1.transpose().dot(vector1));
    }
    else{
        s = index2s[proj_match_point_index] - std::sqrt(vector1.transpose().dot(vector1));
    }
}


    //计算世界坐标系下的x_set， y_set上的点在frenet_path下的坐标s
template <class T>
void world2frenet_path(std::vector<double>& index2s, T& trajectory, msgs::ReferenceLine& rilne, msgs::ReferencePoint& proj_point, int& match_point_index, int n){
    trajectory.frenet_s[0] = CalcSFromIndex2S(index2s, rilne, proj_point, match_point_index);
    Eigen::Vector2d n_r(-std::sin(proj_point.heading), std::cos(proj_point.heading));
    Eigen::Vector2d r_h(trajectory.pose.x - proj_point.x, trajectory.pose.y - proj_point.y);
    trajectory.frenet_l[0] = r_h.transpose().dot(n_r);
    if (n > 0){
        calc_dot_infrenet(trajectory, proj_point);
        if(n > 1){
            calc_dot2_infrenet(trajectory, proj_point);
        }
    }
}

    //计算frenet坐标系下的s_dot， l_dot(frenet_l[1]), dl/ds(frenet_l_dot[0])
template <class T>
void calc_dot_infrenet(T& trajectory, msgs::ReferencePoint& proj_point){
    Eigen::Vector2d v_h(trajectory.velocity.x, trajectory.velocity.y);
    Eigen::Vector2d n_r(-std::sin(proj_point.theta), std::cos(proj_point.theta));
    Eigen::Vector2d t_r(std::cos(proj_point.theta), std::sin(proj_point.theta));
    trajectroy.frenet_l[1] = v_h.transpose().dot(n_r);
    trajectroy.frenet_s[1] = v_h.transpose().dot(t_r/(1 - proj_point.kappa * trajectory.frenet_l[0]));
    if (std::abs(trajectory.frenet_s[1]) < 1e-6){
        trajectory.frenet_l_dot[0] = 0;
    }
    else{
        trajectory.frenet_l_dot[0] = trajectroy.frenet_l[1];
    }
}
    //计算s_dot2, l_dot2, ddl
template <class T>
void calc_dot2_infrenet(T& trajectory, msgs::ReferencePoint& proj_point){
    int i;
    for(i=0; i<trajectory.points.size(); i++){
        if isnan(trajectory.frenet_l[0]){
            break;
        }
        Eigen::Vector2d a_h(trajectory.accel.x, trajectory.accel.y);
        Eigen::Vector2d n_r(-std::sin(proj_point.theta), std::cos(proj_point.theta));
        Eigen::Vector2d t_r(std::cos(proj_point.theta), std::sin(proj_point.theta));
        trajectroy.frenet_l[2] = a_h.transpose().dot(n_r - proj_point * (1 - proj_point.kappa * trajectory.frenet_l[0]);
        trajectroy.frenet_s[2] = (1 / (1 - proj_point.theta * trajectory.frenet_l[0])) * a_h.transpose().dot(t_r) + 2 * proj_point.kappa * trajectory.frenet_l_dot[0] * trajectory.frenet_s[2];
        if (std::abs(trajectory.frenet_s[1]) < 1e-6){
            trajectory.frenet_l_dot[1] = 0;
        }
        else{
            trajectory.frenet_l_dot[1] = (trajectroy.frenet_l[2] - trajectory.frenet_l_dot[0] * trajectory.frenet_s[2]) / std::pow(trajectory.frenet_s[1], 2);
        }
    }
}

template <class T>
double CalcCost(T& pre_point, T& cur_point, int cur_node, auto obs_set){
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

template <class T>
double CalcObsCost(T& square_d, int w_cost_collision){
    double cost = 0;
    for (int i=0; i<square_d.rows(); i++){
        if (square_d(i, 0) < 16 && square_d(i, 0) > 9){
            cost = cost + 1000 / square_d;
        }
        else if (squre_d(i, 0) < 9)
        {
            cost = cost + w_cost_collision;
        }
    }
    return cost
}

void add_density_dp(msgs::Trajectory& trajectory_temp, msgs::TrajectoryPoint& start_point, msgs::Trajectory& trajectory_dp){
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
        auto l = a0 * cal + a1 * s_cur + a2 * cal_pow(s_cur, 2) + a3 * cal_pow(s_cur, 3) + a4 * cal_pow(s_cur, 4) + a5 * cal_pow(s_cur, 5);
        auto dl = a1 * cal + 2 * a2 * s_cur + 3 * a3 * cal_pow(s_cur, 2) + 4 * a4 * cal_pow(s_cur, 3) + 5 * a5 * cal_pow(s_cur, 4);
        auto ddl = 2 * a2 * cal + 6 * a3 * s_cur + 12 * a4 * cal_pow(s_cur, 2) + 20 * a5 * cal_pow(s_cur, 3);
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
        trajectory_dp.points[i].frenet_s[0] = s_cur[i];
        trajectory_dp.points[i].frenet_l[0] = l_temp[i];
        trajectory_dp.points[i].frenet_l[1] = dl_temp[i];
        trajectory_dp.points[i].frenet_l[2] = ddl_temp[i];
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

void add_density_qp(msgs::Trajectory& trajectory_temp,msgs::Trajectory& trajectory_qp){
    int n_init = 20;
    int n = 501;

    double ds = (trajectory_temp.points[-1].frenet_s[0] - trajectory_temp.points[0].frenet_s[0]) / (n-1);
    int index = 1;

    for (int i=0; i<n; i++){
        double x = trajectory_temp.points[-1].frenet_s[0] + i * ds;
        trajectory_qp.points[i].frenet_s[0] = x;
        while (x >= trajectory_temp.points[index].frenet_s[0]){   //while 循环退出的条件是x<qp_path_s(index)，所以x对应的前一个s的编号是index-1 后一个编号是index
            index += 1;
            if (index == n_init){
                break;
            }
        }

        int pre = index - 1;
        int cur = index;
        double delta_s = x - trajectory_temp.points[pre].frenet_s[0];
        double l_pre = trajectory_temp.points[pre].frenet_l[0];
        double dl_pre = trajectory_temp.points[pre].frenet_l[1];
        double ddl_pre = trajectory_temp.points[pre].frenet_l[2];
        double ddl_cur = trajectory_temp.points[cur].frenet_l[2];

        trajectory_qp.points[i].frenet_l[0] = l_pre + dl_pre * delta_s + (1/3) * ddl_pre * std::pow(delta_s, 2) + (1/6)  * ddl_cur * std::pow(delta_s, 2);
        trajectory_qp.points[i].frenet_l[1] = dl_pre + 0.5 * ddl_pre * delta_s + 0.5 * ddl_cur * delta_s;
        trajectory_qp.points[i].frenet_l[2] = ddl_pre + (ddl_cur - ddl_pre) * delta_s / (trajectory_temp.points[cur].frenet_s[0] - trajectory_temp.points[cur].frenet_s[0]);

        index -= 1;
        //因为此时x的后一个编号是index 必有x < qp_path_s(index),在下一个循环中x = x + ds 也未必大于qp_path_s(index)，这样就进入不了while循环，所以index 要回退一位
    }
}