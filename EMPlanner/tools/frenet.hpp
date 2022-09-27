#pragma once
#include <msgs/Object.h>
#include <msgs/ReferencePoint.h>
#include <eigen3/Eigen/Dense>
#include <cmath>

template <class T>
void cartisian2frenet1D(T& object, msgs::ReferencePoint& project_point){
    Eigen::Vector2d pose_vector(object.pose.x - project_point.x, object.pose.y - project_point.y);
    Eigen::Vector2d n_vector(-std::sin(project_point.theta), std::cos(project_point.theta));

    object.frenet_s[0] = project_point.s;
    object.frenet_d_dot[0] = pose_vector.dot(n_vector);
    object.frenet_d[0] = object.frenet_d_dot[0];
}

template <class T>
void cartisian2frenet2D(T& object, msgs::ReferencePoint& project_point){
    Eigen::Vector2d pose_vector(object.pose.x - project_point.x, object.pose.y - project_point.y);
    Eigen::Vector2d velocity_vector(object.velocity.x, object.velocity.y);
    Eigen::Vector2d n_vector(-std::sin(project_point.theta), std::cos(project_point.theta));
    Eigen::Vector2d t_vector(std::cos(project_point.theta), std::sin(project_point.theta));

    object.frenet_s[0] = project_point.s;
    object.frenet_d_dot[0] = pose_vector.dot(n_vector);
    object.frenet_d[0] = object.frenet_d_dot[0];
    object.frenet_s[1] = velocity_vector.dot(t_vector)/(1-project_point.kappa * object.frenet_d[0]);
    object.frenet_d_dot[1] = velocity_vector.dot(n_vector);
    if (std::abs(object.frenet_s[1]) < 1e-6)
        object.frenet_d[1] = 0;
    else
        object.frenet_d[1] = object.frenet_d_dot[1]/object.frenet_s[1];

}

template <class T>
void cartisian2frenet3D(T& object, msgs::ReferencePoint& project_point){
    Eigen::Vector2d pose_vector(object.pose.x - project_point.x, object.pose.y - project_point.y);
    Eigen::Vector2d velocity_vector(object.velocity.x, object.velocity.y);
    Eigen::Vector2d accel_vector(object.accel.x, object.accel.y);
    Eigen::Vector2d n_vector(-std::sin(project_point.theta), std::cos(project_point.theta));
    Eigen::Vector2d t_vector(std::cos(project_point.theta), std::sin(project_point.theta));

    object.frenet_s[0] = project_point.s;
    object.frenet_d_dot[0] = pose_vector.dot(n_vector);
    object.frenet_d[0] = object.frenet_d_dot[0];

    object.frenet_s[1] = velocity_vector.dot(t_vector)/(1-project_point.kappa * object.frenet_d[0]);
    object.frenet_d_dot[1] = velocity_vector.dot(n_vector);
    if (std::abs(object.frenet_s[1]) < 1e-6)
        object.frenet_d[1] = 0;
    else
        object.frenet_d[1] = object.frenet_d_dot[1]/object.frenet_s[1];

    // 默认dkappa=0
    object.frenet_s[2] = (1/(1-project_point.kappa*object.frenet_d[0]))*(accel_vector.dot(t_vector)+2*project_point.kappa * object.frenet_d[1]*object.frenet_s[1]*object.frenet_s[1]);
    object.frenet_d_dot[2] = accel_vector.dot(n_vector) - project_point.kappa*(1-project_point.kappa*object.frenet_d[0]) * std::pow(object.frenet_s[1], 2);
     if (std::abs(object.frenet_s[2]) < 1e-6)
        object.frenet_d[2] = 0;
    else
        object.frenet_d[2] =  (object.frenet_d_dot[2] - object.frenet_d[1]*object.frenet_s[2])/std::pow(object.frenet_s[1], 2);
}

template <class T>
void frenet2cartisian1D(T& object, msgs::ReferencePoint& project_point){
    object.pose.x = project_point.x - sin(project_point.theta) * object.frenet_d[0];
    object.pose.y = project_point.y + cos(project_point.theta) * object.frenet_d[0];
}


