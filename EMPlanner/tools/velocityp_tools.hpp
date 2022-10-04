#pragma once
#include <cmath>
#include <eigen3/Eigen/Dense>
#include "OsqpEigen/OsqpEigen.h"
#include <msgs/Object.h>
#include <msgs/ReferenceLine.h>
#include <msgs/ReferencePoint.h>
#include <tools.hpp>
#include <vector>


template <class T>
void clac_arithmetic(T<double>& vector, double start, double end, double step, int n){
    int length = (end - start) / step + 1;
    cout << length << "some" << endl;
    int j = n - 1;
    vector[0 + j*length] = start;
    for (int i=1; i<length; i++){
        vector[i + j*length] = vector[i-1 + j*length] + step;
    }
}

double CalcObsCost(int s_start, int t_start, int s_end, int t_end, msgs::Velocity_planning_st obs_st, int w_cost_obs){
    double obs_cost = 0;
    int n = 5;
    double dt = (t_end - s_start) / (n - 1);
    double k = (s_end - s_start) / (t_end - t_start);

    for (int i=0; i<n; i++){
        double t = t_start + i * dt;
        double s = s_start + k * i * dt;
        for (int j=0; j<obs_st.s_in.size(); j++){
            Eigen::Vector2d vector1(obs_st.s_in[j] - s, obs_st.t_in[j] - t);
            Eigen::Vector2d vector2(obs_st.s_out[j] - s, obs_st.t_out[j] - t);
            Eigen::Vector2d vector3(obs_st.s_out[j] - obs_st.s_in[j], obs_st.t_out[j] - obs_st.t_in[j]);
            double min_dis = 0;
            double dis1 = std::sqrt(vector1.transpose().dot(vector1));
            double dis2 = std::sqrt(vector2.transpose().dot(vector2));
            double dis3 = std::abs(vector1[0] * vector3[1] - vector1[1] * vector3[0]) / std::sqrt(vector3.transpose().dot(vector3));
            if ((vector1.transpose().dot(vector3) > 0 && vector2.transpose().dot(vector3) > 0) || (vector1.transpose().dot(vector3) <  0 && vector2.transpose().dot(vector3))){
                min_dis = std::min(dis1, dis2);
            }
            else{min_dis = dis3;}
            obs_cost = obs_cost + CalcCollisionCost(w_cost_obs, min_dis);
        }
    }
    return obs_cost
}

double CalcCollisionCost(auto w_cost_obs, auto min_dis){
    double collision_cost;
    if (std::abs(min_dis) < 0.5){
        collision_cost = w_cost_obs;
    }
    else if (std::abs(min_dis) > 0.5 && std::abs(min_dis) < 1.5)
    {
        collision_cost = std::pow(w_cost_obs, ((0.5 - min_dis) + 1));
    }
    else{collision_cost = 0;}
    return collision_cost;
}

double CalcSTCoordinate(int row, int col, auto s_list, auto t_list){  //计算矩阵节点的s t坐标
    int m = s_list.size();
    double s_value = s_list[m - row];
    double t_value = t_list[col];
    return s_value, t_value
}

void add_density_s(msgs::Trajectory& trajectory_temp,msgs::Trajectory& trajectory_qp){
    int t_end = trajectory_temp.points.size();
    auto T = trajectory_temp.points[t_end].t;
    int n = 401;
    auto dt = T / (n - 1);
    for (int i=0; i<n; i++){
        double current_t = (i - 1) * dt;
        for (int j=0; j<t_end-1; j++){
            if (trajectory_temp.points[j].t <= current_t && trajectory_temp.points[j+1] > current_t){
                break;
            }
        }
        auto x = current_t - trajectory_temp.points[j];
        trajectory_qp.points[i].frenet_s[0] = trajectory_temp.points[j].frenet_s[0] + trajectory_temp.points[j].frenet_s[1] * x + (1/3) * trajectory_temp.points[j].frenet_s[2] * std::pow(x, 2) + (1/6) * trajectory_temp.points[j+1].frenet_s[2] * std::pow(x, 2);
        trajectory_qp.points[i].frenet_s[1] = trajectory_temp.points[j].frenet_s[1] 0.5 * trajectory_temp.points[j].frenet_s[2] * x + 0.5 * trajectory_temp.points[j+1].frenet_s[2] * x;
        trajectory_qp.points[i].frenet_s[2] = trajectory_temp.points[j].frenet_s[2] + (trajectory_temp.points[j+1].frenet_s[2] - trajectory_temp.points[j].frenet_s[2]) * x / (trajectory_temp.points[j+1].t - trajectory_temp.points[j].t);
        trajectory_qp.points[i].t = current_t;
    }

}