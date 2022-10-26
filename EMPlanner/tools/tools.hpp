#pragma once
#include <cmath>
#include <eigen3/Eigen/Dense>
#include "OsqpEigen/OsqpEigen.h"
#include <msgs/Object.h>
#include <msgs/ReferencePoint.h>

double cal_distance(double x1, double y1, double x2=0, double y2=0){
    return std::pow(std::pow((x1-x2),2)+std::pow(y1-y2, 2), 0.5);
}

template <class T>
auto cal_pow(T& cal, int n){
    for (auto& i : cal){
        i = std::pow(i, n);
    }
}

template <class T>
double interpolate_linear(T x0, T x1, T y0, T y1, T xi){
    if (xi == x0){
        return y0;
    }
    if (xi == x1){
        return y1
    }
    return y1 * (xi - x0) / (x1 - x0) + y0 * (xi - x1) / (x0 - x1);    
}

template <class T>
void find_match_point(int start, int end, T& object, msgs::Map& map, int& match_point_index, int max_increase_count){
    // 该函数为寻找匹配点的函数，给定寻找的msgs::ReferenceLine类型的路径，以及开始，结束位，
    // 若遍历点到当前距离递增超过50次，则跳出循环
    // 
    // 返回最近点match_point及其match_point_index
    double ex = object.location.pose.x, ey = object.location.pose.y;
    double min_distance = 1e6;
    double pre_distance = 1e6;
    int increase_count = 0;
    if(start<0) start = 0;
    if(end<0) end = 0;
    if(start>map.points.size()) start = map.points.size();
    if(end>map.points.size()) end = map.points.size();
    if(start==end){match_point_index = start;}
    else{
        for(int i=start; i!=end; start<end?i++:i--){
            double distance = cal_distance(ex, ey, map.points[i].rpoint.x, map.points[i].rpoint.y);
            if(distance < min_distance){
                min_distance = distance;
                match_point_index = i;
                increase_count = 0;
            }
            if(pre_distance < distance) increase_count++;
            if(increase_count > max_increase_count) break;
            pre_distance = distance;
        }
    }
    
}

template <class T>
void cal_match_point(int start, int end, T& object, msgs::ReferenceLine& rline, msgs::ReferencePoint& match_point, int& match_point_index, int& max_increase_count){
    // 寻找匹配点, 以rline原点为坐标原点，计算障碍物在参考线上的匹配点match_point（仍为cartesian坐标系下）
    double ex = object.pose.x, ey = object.pose.y;
    double min_distance = 1e6, pre_distance = 1e6;
    int increase_count = 0;
    
    // index2s(msgs::ReferenceLine& rline)

    if(start<0) start = 0;
    if(end<0) end = 0;
    if(start>rline.points.size()) start = rline.points.size();
    if(end>rline.points.size()) end = rline.points.size();
    if(start==end){match_point_index = start;}
    else{
        for(int i=start; i!=end; start<end?i++:i--){
            double distance = cal_distance(ex, ey, rline.points[i].x, rline.points[i].y);
            if(distance < min_distance){
                min_distance = distance;
                match_point_index = i;
                increase_count = 0;
            }
            if(pre_distance < distance) increase_count++;
            if(increase_count > max_increase_count) break;
            pre_distance = distance;
        }
    }
    match_point = rline.points[match_point_index];
}

template <class T>      
void cal_project_point(T& object, msgs::ReferencePoint& match_point, msgs::ReferencePoint& project_point){
    //计算匹配点match_point在frenet坐标轴的投影的直角坐标（proj_x, proj_y, proj_heading, proj_kappa）
    // double ds = object.frenet_s - match_point.s;
    double ds = (object.pose.x - match_point.x)*std::cos(match_point.theta) + (object.pose.y - match_point.y)*std::sin(match_point.theta);
    project_point.x = match_point.x + ds * std::cos(match_point.theta);
    project_point.y = match_point.y + ds * std::sin(match_point.theta);
    project_point.theta = match_point.theta + ds * match_point.kappa;
    project_point.kappa = match_point.kappa;
    project_point.dkappa = match_point.dkappa;
}

template <class T>
Eigen::VectorXd cal_quintic_coef(T&start_point, T&end_point){
    double start_s_1 = start_point.frenet_s[0];
    double start_s_2 = start_s_1*start_s_1;
    double start_s_3 = start_s_2*start_s_1;
    double start_s_4 = start_s_3*start_s_1;
    double start_s_5 = start_s_4*start_s_1;
    double end_s_1 = end_point.frenet_s[0];
    double end_s_2 = end_s_1*end_s_1;
    double end_s_3 = end_s_2*end_s_1;
    double end_s_4 = end_s_3*end_s_1;
    double end_s_5 = end_s_4*end_s_1;

    // a = A^-1 * B
    Eigen::MatrixXd A(6, 6);
    Eigen::VectorXd B(6);
    A << 1,  start_s_1,  start_s_2,   start_s_3,   start_s_4,    start_s_5, 
         0,  1,          2*start_s_1, 3*start_s_2, 4*start_s_3,  5*start_s_4, 
         0,  0,          2,           6*start_s_1, 12*start_s_2, 20*start_s_3, 
         1,  end_s_1,  end_s_2,   end_s_3,   end_s_4,    end_s_5, 
         0,  1,        2*end_s_1, 3*end_s_2, 4*end_s_3,  5*end_s_4, 
         0,  0,        2,         6*end_s_1, 12*end_s_2, 20*end_s_3;
    B << start_point.frenet_l[0], start_point.frenet_l[1], start_point.frenet_l[2], end_point.frenet_l[0], end_point.frenet_l[1], end_point.frenet_l[2];
    return A.inverse() * B;
}

template <class T>
Eigen::VectorXd solve_qp(T H, T f, T A, T b, T Aeq, T beq, T lb, T ub, int n){
    OsqpEigen::Solver solver;
    solver.settings()->setWarmStart(true);
    
    Eigen::SparseMatrix<double> hessian;
    Eigen::VectorXd gradient;
    Eigen::SparseMatrix<double> linearMatrix;
    Eigen::VectorXd lowerBound;
    Eigen::VectorXd upperBound;
    Eigen::MatrixXd Iden = Eigen::MatrixXd::Identity(ub.rows(), ub.rows());
    Eigen::MatrixXd Infinite = Eigen::MatrixXd::Ones(b.rows(), 1);
    Infinite *= -99999;

    int num_Variables = std::max(A.cols(), Aeq.cols());
    int num_Constrains = ub.rows() + Aeq.rows() + A.rows();
    solver.data()->setNumberOfVariables(num_Variables);
    solver.data()->setNumberOfConstraints(num_Constrains);

    hessian.resize(num_Variables, num_Variables);
    gradient.resize(num_Variables);
    linearMatrix.resize(num_Constrains, num_Variables);
    lowerBound.resize(num_Constrains);
    upperBound.resize(num_Constrains);

    hessian.block(0, 0, H.rows(), H.cols()) = H;
    // hessian << H;
    gradient << f;
    linearMatrix.block(0, 0, Iden.rows(), Iden.cols());
    linearMatrix.block(0, Iden.cols(), Aeq.rows(), Aeq.cols());
    linearMatrix.block(0, Iden.cols()+Aeq.cols(), A.rows(), A.cols());
    // linearMatrix << Iden,Aeq,A;
    lowerBound << lb, beq, Infinite;
    upperBound << ub, beq, b;

    if (!solver.data()->setHessianMatrix(hessian)) return false;
    if (!solver.data()->setGradient(gradient)) return false;
    if (!solver.data()->setLinearConstraintsMatrix(linearMatrix)) return false;
    if (!solver.data()->setLowerBound(lowerBound)) return false;
    if (!solver.data()->setUpperBound(upperBound)) return false;
    if (!solver.initSolver()){
        std::cout << "init solver failed" << std::endl;
    }
    if (!solver.solve()){
        std::cout << "qp solve failed" << std::endl;
    }

    Eigen::VectorXd QPSolution;
    QPSolution = solver.getSolution();
    return QPSolution;
}

void CalcProjPoint_of_s(msgs::TrajectoryPoint& trajectory_point, msgs::ReferenceLine& rline, std::vector<double>& index2s, msgs::ReferencePoint& proj_point){
    int match_index = 0;
    double s = trajectory_point.frenet_s[0];
    while (index2s[match_index] < s){
        match_index += 1;
    }
    msgs::ReferencePoint match_point;
    match_point.x = rline.points[match_index].x;
    match_point.y = rline.points[match_index].y;
    match_point.theta = rline.points[match_index].theta;
    match_point.kappa = rline.points[match_index].kappa;

    double ds = s - index2s[match_index];
    proj_point.x = std::cos(match_point.theta) * ds + match_point.x;
    proj_point.y = std::sin(match_point.theta) * ds + match_point.y;
    proj_point.theta = match_point.theta + ds * match_point.kappa;
    proj_point.kappa = match_point.kappa;
}

void path_frenet2cartesian(msgs::Trajectory& trajectory, msgs::ReferenceLine& rline, std::vector<double>& index2s){
    msgs::ReferencePoint proj_point;
    for (auto trajectory_point : trajectory.points){
        CalcProjPoint_of_s(trajectory_point, rline, index2s, proj_point);
        trajectory_point.pose.x = proj_point.x + trajectory_point.frenet_l[0] * -std::sin(proj_point.theta);
        trajectory_point.pose.y = proj_point.y + trajectory_point.frenet_l[0] * std::cos(proj_point.theta);
        trajectory_point.pose.theta = proj_point.theta + std::atan(trajectory_point.frenet_l[1] / (1 - proj_point.kappa * trajectory_point.frenet_l[0]));
        trajectory_point.kappa = ((trajectory_point.frenet_l[2] + proj_point.theta * trajectory_point.frenet_l[1] * std::tan(trajectory_point.pose.theta - proj_point.theta)) * std::pow(std::cos(trajectory_point.pose.theta - \
            proj_point.theta), 2) / (1 - proj_point.kappa * trajectory_point.frenet_l[0]) + proj_point.kappa) * std::cos(trajectory_point.pose.theta - proj_point.theta) / (1 - proj_point.kappa * trajectory_point.frenet_l[0]);
    }
}


void path_merge_velocity(msgs::TrajectoryPoint& start_point, msgs::Trajectory& path, msgs::Trajectory& velocity, msgs::Trajectory& trajectory_final){
    int index = path.points.size() - 1;
    int n = 401;
    auto current_time = start_point.t;

    // std::vector<auto> temp;
    // std::vector<auto>::const_iterator first = path.points.begin();
    // std::vector<auto>::const_iterator second = path.points.begin()+index;
    // temp.assign(first, second);
    // std::vector<int> temp{&path.points[i], &path.points[i]+index};

    for (int i=0; i<n; i++){
        trajectory_final.points[i].pose.x = interpolate_linear(path.points[0].frenet_s[0], path.points[index].frenet_s[0], path.points[0].pose.x, path.points[index].pose.x, velocity.points[i].frenet_s[0]);
        trajectory_final.points[i].pose.y = interpolate_linear(path.points[0].frenet_s[0], path.points[index].frenet_s[0], path.points[0].pose.y, path.points[index].pose.y, velocity.points[i].frenet_s[0]);
        trajectory_final.points[i].pose.theta = interpolate_linear(path.points[0].frenet_s[0], path.points[index].frenet_s[0], path.points[0].pose.theta, path.points[index].pose.theta, velocity.points[i].frenet_s[0]);
        trajectory_final.points[i].kappa = interpolate_linear(path.points[0].frenet_s[0], path.points[index].frenet_s[0], path.points[0].kappa path.points[index].kappa, velocity.points[i].frenet_s[0]);
        trajectory_final.points[i].t = velocity.points[i].t + current_time;
        trajectory_final.points[i].velocity.x = velocity.points[i].frenet_s[1];
        trajectory_final.points[i].accel.x = velocity.points[i].frenet_s[2];
    }
}

void stitch_trajectory(msgs::Trajectory& init, msgs::Trajectory& stitch, msgs::Trajectory& final){     //此处的stitch是路径规划计算起点的stitch
    final.points.insert(final.points.begin(), init.points.begin(), init.points.end());
    final.points.insert(final.points.end(), stitch.points.begin(), stitch.points.end());
}