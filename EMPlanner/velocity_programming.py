from re import L
import numpy as np
import math
from cmath import atan, isnan, nan, sqrt

REFERRENCE_VELOCITY = 15

class vp:
    def __init__(self) -> None:
        pass

        #计算trajectory的s与x, y的对应关系，可以看作是trajecotry index2s
    def trajectory_s2xy(self, trajectory_x, trajectory_y):
        n = len(trajectory_x)
        path_index2s = np.zeros((n, 1))
        sum = 0
        
        for i in range(1, len(trajectory_x)):
            if isnan (trajectory_x[i]):
                break
            
            sum += sqrt((trajectory_x[i] - trajectory_x[i-1])**2 + (trajectory_y[i] - trajectory_y[i-1])**2)
            path_index2s[i] = sum

        if i == n:
            path_s_end = path_index2s[-1]
        else:
            path_s_end = path_index2s[i-1]
        
        return path_s_end, path_index2s

        #计算速度规划的起始条件
    def programming_init(self, plan_start_vx, plan_start_vy, plan_start_ax, plan_start_ay, plan_start_heading):
        tor = [math.cos(plan_start_heading), math.sin(plan_start_heading)]
        
        v_t = np.dot(np.matrix(tor).T, [plan_start_vx, plan_start_vy])  #计算向量v，a 在切向的投影
        a_t = np.dot(np.matrix(tor).T, [plan_start_ax, plan_start_ay])
        plan_start_s_dot, plan_start_s_dot2 = v_t, a_t
        return plan_start_s_dot, plan_start_s_dot2

        #生成st图，使用斜直线模型
    def generate_st_graph(self, pre_trajectory_s, pre_trajectory_t, dynamic_obs_s_set, dynamic_obs_l_set, dynamic_obs_s_dot_set, dynamic_obs_l_dot_set):
        n = len(dynamic_obs_s_set)
         
        obs_st_s_in_set, obs_st_s_out_set, obs_st_t_in_set, obs_st_t_out_set, obs_st_l_in_set, obs_st_l_out_set = np.ones((n, 1)) * nan, np.ones((n, 1)) * nan, np.ones((n, 1)) * nan,  np.ones((n, 1)) * nan, np.ones((n, 1)) * nan, np.ones((n, 1)) * nan
        for i in range(len(dynamic_obs_s_set)):
            if isnan(dynamic_obs_s_set[i]):
                break
            if abs(dynamic_obs_l_dot_set[i]) < 0.3:
                if abs(dynamic_obs_l_set[i]) > 2:   #横向距离较远，且横向速度缓慢，可忽略
                    continue
                else:   #移动缓慢且距离较近
                    #TODO需要做虚拟障碍物，感知模块加入跟踪逻辑，判断两帧之间的障碍物是否是同一个，给到速度规划模块后，先给出虚拟障碍物觉，下一帧拿到虚拟障碍物的信息，规划绕过去的路径和速度
                    #本算法欠缺的：障碍物结构体，包含了坐标、速度、还要包含决策标记（是否为虚拟障碍物，左绕还是右绕，避让还是超车
                    continue
            t_zero = -dynamic_obs_l_set[i] / dynamic_obs_l_dot_set[i] #动态障碍物的坐标l到与原点0所需要的时间
            t_boundary1 = 2 / dynamic_obs_l_dot_set[i] + t_zero #加上缓冲时间
            t_boundary2 = -2 / dynamic_obs_l_dot_set[i] + t_zero
            t_max = max(t_boundary1, t_boundary2)
            t_min = min(t_boundary1, t_boundary2)
            
            if t_max < 1 or t_min > 8:  #忽略太远或者太近（碰瓷的）的
                #遇上碰瓷的，即使做出很大加速度也执行不了， 碰瓷的障碍物也需要做虚拟障碍物，和路径一起解决
                continue

            if t_min < 0 and t_max > 0: #感知看到时，障碍物已经在±2的内部了
                obs_st_s_in_set[i] = dynamic_obs_s_set[i]
                obs_st_s_out_set[i] = dynamic_obs_s_set[i] + dynamic_obs_s_dot_set[i] * t_max
                obs_st_t_in_set[i] = 0
                obs_st_t_out_set[i] = t_max
            else:   #正常障碍物
                obs_st_s_in_set[i] = dynamic_obs_s_set[i] + dynamic_obs_s_dot_set[i] * t_min
                obs_st_s_out_set[i] = dynamic_obs_s_set[i] + dynamic_obs_s_dot_set[i] *t_max
                obs_st_t_in_set[i] = t_min
                obs_st_t_out_set[i] = t_max


        return obs_st_s_in_set, obs_st_s_out_set, obs_st_t_in_set, obs_st_t_out_set

        #速度的动态规划
    def speed_dp(self, obs_st_s_in_set, obs_st_s_out_set, obs_st_t_in_set, obs_st_t_out_set, reference_speed_unlimit, w_cost_ref_speed, w_cost_accel, w_cost_obs, plan_start_s_dot):
        #时间从0到8开始规划，最多8秒 ； s的范围从0开始到路径规划的path的总长度 ； 为减少算力，采用非均匀采样，越小的越密，越大的越稀疏
        s_list = [np.arange(start=0, stop=4.5, step=0.5), np.arange(start=5.5, stop=14.5, step=1), np.arange(start=16, stop=29.5, step=1.5), np.arange(start=32, stop=54.5, step=2.5)]
        s_list = np.reshape(s_list, (len(s_list, s_list[0])))
        t_list = np.arange(start=0.5, stop=8, step=0.5)
            #s采样40个点，t采样16个点
        dp_st_cost = np.ones((40, 16)) * 99999  #代价矩阵
        dp_st_node = np.ones((40, 16))
        dp_st_s_dot = np.zeros((40, 16))        #表示从起点到i，j点的最优路径的末速度

        for i in range(len(s_list)):    #计算从dp起点到第一列的cost
            dp_st_cost[i][1] = self.CalcDPCost(0, 0, i, 1, obs_st_s_in_set, obs_st_s_out_set, obs_st_t_in_set, obs_st_t_out_set, w_cost_ref_speed, reference_speed, w_cost_accel, w_cost_obs, plan_start_s_dot, s_list, t_list, dp_st_s_dot)
            s_end, t_end = self.CalcSTCoordinate(i, 1, s_list, t_list)
            dp_st_s_dot[i][1] = s_end / t_end
        
        for i in range(1, len(t_list)):
            for j in range(len(s_list)):
                cur_row = j
                cur_col = i
                for k in range(len(s_list)):
                    pre_row = k
                    pre_col = i - 1
                    cost_temp = self.CalcDPCost(self, pre_row, pre_col, cur_row, cur_col, obs_st_s_in_set, obs_st_s_out_set, obs_st_t_in_set, obs_st_t_out_set, w_cost_ref_speed, reference_speed, w_cost_accel, w_cost_obs, plan_start_s_dot, s_list, t_list, dp_st_s_dot)

                    if cost_temp + dp_st_cost[pre_row][pre_col] < dp_cost[cur_row][cur_col]:
                        dp_cost[cur_row][cur_col] = cost_temp + dp_st_cost[pre_row, pre_col]
                        #计算最优的s_dot
                        s_start, t_start = self.CalcSTCoordinate(pre_row, pre_col, s_list, t_list)
                        s_end, t_end = self.CalcSTCoordinate(cur_row, cur_col, s_list, t_list)
                        dp_st_s_dot[cur_row][cur_col] = (s_end - s_start) / (t_end - t_start)
                        #将最短路径的前一个节点的行号记录下来
                        dp_st_node[cur_row][cur_col] = pre_row
        
        dp_speed_s = np.ones((len(t_list, 1))) * nan
        dp_speed_t = dp_speed_s
        #找到dp_node_cost 上边界和右边界代价最小的节点
        min_cost, min_row, min_col = 99999, 99999, 99999
        for i in range(len(s_list)):    #遍历右边界
            if dp_st_cost[i][len(t_list)] < min_cost:
                min_cost = dp_st_cost[i][len(t_list)]
                min_row = i
                min_col = len(t_list)
            
        for j in range(len(t_list)):    #遍历上边界
            if dp_st_cost[1][j]  <= min_cost:
                min_cost = dp_st_cost
                min_row = 1
                min_col = j
            
        #这里要注意，虽然我们在t上每0.5s 采样一个时间，采了16个点，但是min_col未必等于16， 也就是说动态规划的最优解的时间未必是8秒
        dp_speed_s[min_col], dp_speed_t[min_col] = self.CalcSTCoordinate(min_row, min_col, s_list, t_list)  #先计算终点
        while min_col != 1:     #回溯
            pre_row = dp_st_node[min_row][min_col]
            pre_col = min_col - 1
            dp_speed_s[pre_col], dp_speed_t[pre_col] = self.CalcSTCoordinate(pre_row, pre_col, s_list, t_list)
            min_row, min_col = pre_row, pre_col
        
        return dp_speed_s, dp_speed_t


        #计算连接两个节点之间的边的代价
    def CalcDPCost(self, row_start, col_start, row_end, col_end, obs_st_s_in_set, obs_st_s_out_set, obs_st_t_in_set, obs_st_t_out_set, w_cost_ref_speed, reference_speed, w_cost_accel, w_cost_obs, plan_start_s_dot, s_list, t_list, dp_st_s_dot):
        #输入：边的起、终点的行列号：row_start, col_start, row_end, col_end ； 障碍物st信息：obs_st_s_in_set, obs_st_s_out_set, obs_st_t_in_set, obs_st_t_out_set ；
        # 代价权重w_cost_ref_speed（推荐速度代价, w_cost_accel(加速度代价权重), w_cost_obs(障碍物代价) ； 推荐速度reference_speed ； 拼接起点的速度plan_start_s_dot ； 采样距离s、t list ； 上节点最优路径的末速度dp_st_s_dot 用于计算加速度

        s_end, t_end = self.CalcSTCoordinate(row_end, col_end, s_list, t_list)
        if row_start == 0:  #边的起点是dp的起点
            s_start = 0
            t_start = 0
            s_dot_start = plan_start_s_dot
        else:
            s_start, t_start = self.CalcSTCoordinate(row_start, col_start, s_list, t_list)
            s_dot_start = dp_st_s_dot[row_start, col_start]
        
        cur_s_dot = (s_end - s_start) / (t_end - t_start)
        cur_s_dot2 = (cur_s_dot - s_dot_start) / (t_end - t_start)
        cost_ref_speed = w_cost_ref_speed * (cur_s_dot - reference_speed)**2
        if cur_s_dot2 < 4 and cur_s_dot2 > -6:
            cost_accel = w_cost_accel * cur_s_dot2**2
        else:   #超过了车辆动力学上下限
            cost_accel = 9999 * w_cost_accel * cur_s_dot2**2
        cost_obs = self.CalcObsCost(s_start, t_start, s_end, t_end, obs_st_s_in_set, obs_st_s_out_set, obs_st_t_in_set, obs_st_t_out_set, w_cost_obs)
        cost = cost_obs + cost_accel + cost_ref_speed

        return cost

    def CalcObsCost(self, s_start, t_start, s_end, t_end, obs_st_s_in_set, obs_st_s_out_set, obs_st_t_in_set, obs_st_t_out_set, w_cost_obs):
        #输入：边的起点终点：s_start, t_start, s_end, t_end ； 障碍物信息：obs_st_s_in_set, obs_st_s_out_set, obs_st_t_in_set, obs_st_t_out_set ； 障碍物代价：w_cost_obs
        obs_cost = 0
        n =5    #采样点个数
        dt = (t_end - s_start) / (n - 1)    #采样时间间隔
        k = (s_end - s_start) / (t_end - t_start)   #边的斜率
        for i in range(n):
            t = t_start + i * dt    #计算采样点的坐标
            s = s_start + k * i * dt
            for j in range(len(obs_st_s_in_set[j])):    #遍历所有障碍物
                if isnan(obs_st_s_in_set[j]):
                    continue
                #计算点到st线段的距离
                vector1 = [obs_st_s_in_set[j], obs_st_t_in_set[j]] - [s, t]
                vector2 = [obs_st_s_out_set[j], obs_st_t_out_set[j]] - [s, t]
                vector3 = vector2 - vector1
                min_dis = 0
                dis1 = sqrt(np.dot(vector1.T, vector1))
                dis2 = sqrt(np.dot(vector2.T, vector2))
                dis3 = abs(np.dot(vector1[0], vector3[1]) - np.dot(vector1[1], vector3[0])) / sqrt(np.dot(vector3.T, vector3))
                if (np.dot(vector1.T, vector3) > 0 and np.dot(vector2.T, vector3) > 0) or (np.dot(vector1.T, vector3) < 0 and np.dot(vector2.T,vector3) < 0):
                    min_dis = min(dis1)
                else:
                    min_dis = dis3
                obs_cost = obs_cost +  self.CalcCollisionCost(w_cost_obs, min_dis)

        return obs_cost

    def CalcCollisionCost(self, w_cost_obs, min_dis):
        if abs(min_dis) < 0.5:
            collision_cost = w_cost_obs
        elif abs(min_dis) > 0.5 and abs(min_dis) < 1.5:
            collision_cost = w_cost_obs**((0.5 - min_dis) + 1)
        else:
            collision_cost = 0

        return collision_cost

        #计算矩阵节点的s t坐标
    def CalcSTCoordinate(self, row, col, s_list, t_list):
        #输入：row col 节点在矩阵的行号和列号， s_list t_list 采样间隔表
        #输出：s_value, t_value节点的s t坐标
        #矩阵的（1,1）代表的是最左上角的元素 s最大， t最小

        m = len(s_list)
        s_value = s_list[m - row]
        t_value = t_list[col]
        return s_value, t_value