#pragma once
#include <cmath>
#include <eigen3/Eigen/Dense>

double cal_distance(double x1, double y1, double x2=0, double y2=0){
    return std::pow(std::pow((x1-x2),2)+std::pow(y1-y2, 2), 0.5);
}

//     //计算index与s的转换，index2s（i）表示编号i对应的弧长s
// template <class T>
// void index2s(msgs::ReferenceLine& rline, msgs::ReferencePoint& origin_point, int& origin_match_point_index){
//     for (int i=1; i<=rline.points.size(); i++){
//         rline.points[i].s = std::sqrt(std::pow(rline.points[i].x - rline.points[i-1].x, 2) + std::pow(rline[i].points.y - rline[i-1].points.y, 2)) + rline.points[i].s[i-1];
//     }
    
//     s0 = CalcSfromIndex2s(index2s, rline, origin_point, origin_match_point_index);
//     for (int i=0; i<=180; i++){
//         trajectory.points[i] = index2s[i] - s0;
//     }
// }

//     //计算点proj_x， proj_y对应的弧长，并判断投影点在匹配点的前面还是后面
// template <class T>
// void CalcSfromIndex2s(int& index2s[181], msgs::ReferenceLine& rline, msgs::ReferencePoint& project_point, int& proj_match_point_index){
//     Eigen::Vector2d vector1(project_point.x - rline.points[proj_match_point_index].x, project_point.y - rline.points[proj_match_point_index].y);
//     if (proj_match_point_index < rline.points.size()){
//         Eigen::Vector2d vector2(rline.points[proj_match_point_index + 1].x - rline.points[proj_match_point_index].x, rline.points[proj_match_point_index + 1].y - rline.points[proj_match_point_index].y);
//     }
//     else{
//         Eigen::Vector2d vector2(rline.points[proj_match_point_index].x - rline.points[proj_match_point_index - 1].x, rline.points[proj_match_point_index].y - rline.points[proj_match_point_index - 1].y);
//     }
    
//     if (vector1.transpose().dot(vector2) > 0){      //投影点在匹配点前面
//         s = index2s[proj_match_point_index] + std::sqrt(vector1.transpose().dot(vector1));
//     }
//     else{
//         s = index2s[proj_match_point_index] - std::sqrt(vector1.transpose().dot(vector1));
//     }
// }

template <class T>
void find_match_point(int start, int end, T& object, msgs::ReferenceLine& rline, msgs::ReferencePoint& match_point, int& match_point_index, int& max_increase_count){
    // 寻找匹配点, 以rline原点为坐标原点，计算障碍物在参考线上的匹配点match_point（仍为cartesian坐标系下）
    double ex = object.location.pose.x, ey = object.location.pose.y;
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
    match_point = rilne.points[match_point_index];
}

template <class T>      
void cal_project_point(T& object, msgs::ReferencePoint& match_point, msgs::ReferencePoint& project_point, int& match_point_index, int[]& index2s){
    //计算匹配点match_point在frenet坐标轴的投影的直角坐标（proj_x, proj_y, proj_heading, proj_kappa）
    double ds = object.location.frenet_s - match_point.s;      
    // double ds = (object.pose.x - match_point.x)*std::cos(match_point.theta) + (object.pose.y - match_point.y)*std::sin(match_point.theta);
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
    B << start_point.frenet_d[0], start_point.frenet_d[1], start_point.frenet_d[2], end_point.frenet_d[0], end_point.frenet_d[1], end_point.frenet_d[2];
    return A.inverse() * B;
}

template <class T>
void CalcObsCost(T& square_d, int w_cost_collision){
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