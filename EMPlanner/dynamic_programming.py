from cmath import atan, isnan, sqrt
import numpy as np
import math


W_COST_COLLISION = 999999   #权重参数
W_COST_SMOOTH_DL = 300
W_COST_SMOOTH_DDL = 2000
W_COST_SMOOTH_DDDL = 10000
W_COST_REF = 20

ROW = 9     #九行六列
COL = 6
SAMPLE_S = 10   #每十米采样
SAMPLE_L = 1    #每一米采样

class Dynamic_Planning:
    def __init__(self) -> None:
        self.is_first_run = 0 
        pass

    def transposition(self, matrix):
        return [[matrix[i][j] for i in range(len(matrix))] for j in range(len(matrix[0]))]

    #计算规划起点以及拼接轨迹的信息
    def start_point_and_stitch_trajectory(self, pre_trajectory_x, pre_trajectory_y, pre_trajectory_heading, pre_trajectory_kappa, pre_trajectory_velocity, pre_trajectory_accel, pre_trajectory_time, current_time, 
                                            host_x, host_y, host_heading_xy, host_kappa, host_vx, host_vy, host_heading, host_ax, host_ay): 
        #输入为：上一周期的规划结果信息pre_trajectory x, y, heading, kappa, velocity, accel, time ; 本周期的绝对时间current_time ; 定位信息host x, y, heading , kappa, vx, vy, ax, ay 
        #输出为：规划起点信息plan_start x, y, heading, kappa, vx, vy, ax, ay(世界坐标系); 待拼接的轨迹信息， 一般在上一周期的轨迹中的curent_time +100 ms, 往前取20个点。
        #规划完成后， 本周期的规划结果和stitch_trajectory一起拼好发给控制
        
        #待拼接轨迹初始化（输出初始化）
        stitch_x, stitch_y = np.zeors((20, 1)), np.zeors((20, 1))
        stitch_heading, stitch_kappa, stitch_speed, stitch_accel = np.zeros((20, 1)), np.zeors((20, 1)), np.zeors((20, 1)), np.zeors((20, 1))
        stitch_time = np.ones((20, 1))  #时间为负代表无拼接轨迹
        
        if self.is_first_run:   #代码首次运行，无上周期轨迹
            self.is_first_run = 0
            #车刚启动，无速度，没必要动力学递推， 直接用定位信息
            plan_start_x, plan_start_y = host_x, host_y
            plan_start_heading, plan_start_kappa, plan_start_vx, plan_start_vy, plan_start_ax, plan_start_ay = host_heading, 0, 0, 0, 0, 0
            plan_start_time  = current_time + 0.1 #规划起点时间应该是当前时间加100ms

        else:
            x_cur, y_cur = host_x, host_y
            heading_cur, kappa_cur = host_heading_xy, 0

            #vx, vy, ax, ay要转换成世界坐标系
            vx_cur = host_vx * math.cos(heading_cur) - host_vy * math.sin(heading_cur)
            vy_cur = host_vx * math.sin(heading_cur) + host_vy * math.cos(heading_cur)
            ax_cur = host_ax * math.cos(heading_cur) - host_ay * math.sin(heading_cur)
            ay_cur = host_ax * math.sin(heading_cur) + host_ay * math.cos(heading_cur)

            dt = 0.1 #规划周期
            
            #上一周期轨迹的本周期车辆理论位置
            for cur_t in range(len(pre_trajectory_time)):
                if pre_trajectory_time[cur_t] <= current_time and pre_trajectory_time[cur_t+1] > current_time:
                    break
            pre_x_desire, pre_y_desire, pre_heading_desire = pre_trajectory_x[cur_t], pre_trajectory_y[cur_t], pre_trajectory_heading[cur_t]
            #计算横纵向误差
            tor = [math.cos(pre_heading_desire), math.sin(pre_heading_desire)]   #切向量
            nor = [-math.sin(pre_heading_desire), math.cos(pre_heading_desire)]     #法向量
            d_err = [host_x, host_y] - [pre_x_desire, pre_y_desire]             #误差向量
            d_err_t = self.transposition(d_err)  #转置
            lon_err, lat_err = abs(np.dot(d_err_t, tor)), abs(np.dot(d_err_t, nor))  #纵、横向误差
            
            if  lon_err > 2.5 or lat_err > 0.5: #若误差过大，动力学外推
                plan_start_x = x_cur + vx_cur * dt + 0.5 * ax_cur * dt * dt
                plan_start_y = y_cur + vy_cur * dt + 0.5 * ay_cur * dt * dt
                plan_start_vx = vx_cur + ax_cur * dt
                plan_start_vy = vy_cur + ay_cur *dt
                plan_start_heading = math.atan2(plan_start_vy, plan_start_vx)
                plan_start_ax, plan_start_ay, plan_start_kappa = ax_cur, ay_cur, kappa_cur
                plan_start_time = current_time + 0.1
            else:                   #控制正常，拼接
                for j in range(cur_t, len(pre_trajectory_time) - 1):
                    if pre_trajectory_time[j] <= current_time + 0.1 and pre_trajectory_time[j+1] <= current_time + 0.1:
                        break
                plan_start_x = pre_trajectory_x[j]  #算法对轨迹点的密度有一定要求
                plan_start_y = pre_trajectory_y[j]
                plan_start_heading = pre_trajectory_heading[j]
                plan_start_kappa = pre_trajectory_kappa[j]
                plan_start_vx = pre_trajectory_velocity[j] * math.cos(plan_start_heading)   #这里的pre_trajectory的速度、加速度是指轨迹的切向速度和切向加速度
                plan_start_vy = pre_trajectory_velocity[j] * math.sin(plan_start_heading)
                tor = [math.cos(plan_start_heading), math.sin(plan_start_heading)]
                nor = [-math.sin(plan_start_heading), math.cos(plan_start_heading)]
                a_tor = pre_trajectory_accel[j] * tor
                a_nor = pre_trajectory_accel[j]**2 * plan_start_kappa * nor
                plan_start_ax = a_tor[0] + a_nor[0]
                plan_start_ay = a_tor[1] + a_nor[1]
                plan_start_time = pre_trajectory_time[j]
                
                #计算拼接轨迹， j代表current_time + 0.1的时间点在pre_trajectory_time
                #拼接从j开始，往前拼接20个，也就是pre_trajectory[j-1]、[j-2]...[j-19]，因此要考虑j是否小于二十
                #pre_trajectory_x,_y 是规划的起点，那么轨迹拼接的时候不能包含进去该点，因为在规划模块中已经包含，拼接时再包含会重复
                j = j -1    #移除规划起点
                if j >= 20:
                    stitch_x = pre_trajectory_x[j-19: j+1]
                    stitch_y = pre_trajectory_y[j-19: j+1]
                    stitch_heading = pre_trajectory_heading[j-19: j+1]
                    stitch_kappa = pre_trajectory_kappa[j-19: j+1]
                    stitch_speed = pre_trajectory_velocity[j-19: j+1]
                    stitch_accel = pre_trajectory_accel[j-19: j+1]
                    stitch_time = pre_trajectory_time[j-19: j+1]
                else:
                    stitch_x[20-j+1, 21] = pre_trajectory_x[1: j+1]
                    stitch_y[20-j+1, 21] = pre_trajectory_y[1: j+1]
                    stitch_heading[20-j+1, 21] = pre_trajectory_heading[1: j+1]
                    stitch_kappa[20-j+1, 21] = pre_trajectory_kappa[1: j+1]
                    stitch_speed[20-j+1, 21] = pre_trajectory_velocity[1: j+1]
                    stitch_accel[20-j+1, 21] = pre_trajectory_accel[1: j+1]
                    stitch_time[20-j+1, 21] = pre_trajectory_time[1: j+1]

    #计算index与s的转换，index2s（i）表示编号i对应的弧长s
    def index2s(self, path_x, path_y, origin_x, origin_y, origin_match_point_index):
        n = 181
        index2s = np.zeros((n,1))
        
        for i in range(1, 181):#计算以轨迹起点为坐标原点的index2s的转换关系
            index2s[i] = math.sqrt((path_x[i] - path_x[i-1])**2 + (path_y[i] - path_y[i-1])**2) + index2s[i-1]
        
        s0 = self.CalcSFromIndex2S(index2s, path_x, path_y, origin_x, origin_y, origin_match_point_index)
        index2s = index2s - np.ones((n, 1))*s0

    #给定index2s的映射关系后，计算点（proj_x， proj_y）所对应的弧长，并判断投影点在匹配点的前面还是后面（向量内积
    def CalcSFromIndex2S(self, index2s, path_x, path_y, proj_x, proj_y, proj_match_point_index):
        vector_1 = [proj_x, proj_y] - [path_x[proj_match_point_index], path_y[proj_match_point_index]]
        vector_1_t = self.transposition(vector_1)
        if proj_match_point_index < len(path_x):
            vector_2 = [path_x[proj_match_point_index + 1], path_y[proj_match_point_index + 1]] - [path_x[proj_match_point_index], path_y[proj_match_point_index]]
        else:
            vector_2 = [path_x[proj_match_point_index], path_y[proj_match_point_index]] - [path_x[proj_match_point_index - 1], path_y[proj_match_point_index - 1]]

        if np.dot(vector_1_t, vector_2) > 0:   #t投影点在匹配点前面s
            s = index2s[proj_match_point_index] + sqrt(np.dot(vector_1_t, vector_1))
        else:
            s = index2s[proj_match_point_index] - sqrt(np.dot(vector_1_t, vector_1))

    #计算世界坐标系下的x_set， y_set上的点在frenet_path下的坐标s
    def world2frenet_path(self, x_set, y_set, frenet_path_x, frenet_path_y, proj_x_set, proj_y_set, proj_heading_set, proj_match_point_index_set, index2s):
        #输入为： x_set, y_set,待转换的点， frenet_path_x, frenet_path_y, frenet坐标点， proj_x_set, proj_y_set, porj_heading_set, proj_match_point_index_set, index2s 待坐标转换的点的投影点
        n = 128 #未知需要坐标转换点的个数，直接用最大值

        s_set = np.ones((n, 1))
        l_set = np.ones((n, 1))
        for i in range(len(x_set)):
            if isnan(x_set[i]):
                break
            s_set[i] = self.CalcSFromIndex2S(index2s, frenet_path_x, frenet_path_y, proj_x_set[i], proj_y_set[i], proj_match_point_index_set[i])
            n_r = [-math.sin(proj_heading_set[i]), math.cos(proj_heading_set[i])]
            r_h = [x_set[i], y_set[i]]
            r_r = [proj_x_set[i], proj_y_set[i]]
            l_set[i] = np.dot(self.transposition(r_h - r_r), n_r
)
    #计算frenet坐标系下的s_dot， l_dot, dl/ds
    def cal_dot(self, l_set, vx_set, vy_set, proj_heading_set, proj_kappa_set):
        n = 128

        s_dot_set = np.ones((n, 1))
        l_dot_set = np.ones((n, 1))
        dl_set = np.ones((n ,1))
 
        for i in range(len(l_set)):
            if isnan(l_set[i]):
                break
            v_h = [vx_set[i], vy_set[i]]
            n_r = [-math.sin(proj_heading_set[i]), math.cos(proj_heading_set[i])]
            t_r = [math.cos(proj_heading_set[i]), math.sin(proj_heading_set[i])]
            l_dot_set[i] = np.dot(self.transposition(v_h), n_r)
            s_dot_set[i] = np.dot(self.transposition(v_h), t_r /(1 - proj_kappa_set[i] * l_set[i]))
            if abs(s_dot_set[i]) < 1e-6:
                dl_set[i] = 0
            else:
                dl_set[i] = l_dot_set [i]

    #计算s_dot2, l_dot2, ddl
    def cal_dot2(self, ax_set, ay_set, proj_heading_set, proj_kappa_set, l_set, s_dot_set, dl_set):
        n = 128

        s_dot2_set = np.ones((n, 1))
        l_dot2_set = np.ones((n, 1))
        ddl_set = np.ones((n, 1))

        for i in range(len(l_set)):
            if isnan(l_set[i]):
                break
            a_h = [ax_set[i], ay_set[i]]
            n_r = [-math.sin(proj_heading_set[i]), math.cos(proj_heading_set[i])]
            t_r = [math.cos(proj_heading_set[i]), math.sin(proj_heading_set[i])]
            l_dot2_set[i] = np.dot(self.transposition(a_h), n_r - proj_kappa_set[i] * (1 - proj_kappa_set[i] * l_set))
            s_dot2_set[i] = (1 / (1 - proj_heading_set[i] * l_set[i])) * np.dot((self.transposition(a_h), t_r) + 2 * proj_kappa_set[i] * dl_set[i] * s_dot2_set[i])
            if abs(s_dot2_set[i]) < 1e-6:
                ddl_set[i] = 0
            else:
                ddl_set[i] = (l_dot2_set[i] - dl_set[i] * s_dot2_set[i]) / (s_dot_set[i]**2)

    #过滤静态障碍， 仅单车道
    def deal_with_obstacle(self, host_x, host_y, host_heading_xy, obs_x_set_gcs, obs_y_set_gcs, obs_velocity_set_gcs, obs_heading_set_gcs):
        n = 32 #传感器最多可识别的障碍物
        count = 1

        obs_x_set_final = np.ones((n, 1))
        obs_y_set_final = np.ones((n, 1))
        obs_velocity_set_final = np.ones((n, 1))
        obs_heading_set_final = np.ones((n, 1))

        for i in range(len(obs_x_set_gcs)):
            if isnan(obs_x_set_gcs[i]):
                break
            tor = [math.cos(host_heading_xy), math.sin(host_heading_xy)]    #自车heading的方向向量与法向量
            nor = [-math.sin(host_heading_xy), math.cos(host_heading_xy)]   
            vector_obs = [obs_x_set_gcs[i], obs_y_set_gcs[i]] - [host_x, host_y]    #障碍物与自车的距离向量
            lon_distance = np.dot(self.transposition(vector_obs), tor) #障碍物与自车距离
            lat_distance = np.dot(self.transposition(vector_obs), nor)
            if lon_distance < 60 and lon_distance > -10 and lat_distance > -10 and lat_distance < 10:
            # if lon_distance < 60 and lon_distance > -10 and lat_distance > -25 and lat_distance < 25: #加入动态障碍物时
                obs_x_set_final[count] = obs_x_set_gcs[i]
                obs_y_set_final[count] = obs_y_set_gcs[i]
                obs_heading_set_final[count] = obs_heading_set_gcs[i]
                obs_velocity_set_final[count] = obs_velocity_set_gcs[i]
                count += 1

    #筛选分类静态障碍物与动态障碍物
    def distinguish_obstacle(self, obs_x_set_gcs, obs_y_set_gcs, obs_heading_set_gcs, obs_velocity_set_gcs):
        n = 32
        static_obs_x_set = np.ones((n, 1))
        static_obs_y_set = np.ones((n, 1))
        dynamic_obs_x_set = np.ones((n, 1))
        dynamic_obs_y_set = np.ones((n, 1))
        dynamic_obs_vx_set = np.ones((n, 1))
        dynamic_obs_vy_set = np.ones((n, 1))
        count_static = 1
        count_dynamic = 1

        for i in range(len(obs_x_set_gcs)):
            if abs(obs_velocity_set_gcs[i]) < 0.1:
                static_obs_x_set[count_static] = obs_x_set_gcs[i]
                static_obs_y_set[count_static] = obs_y_set_gcs[i]
                count_static += 1
            else:
                dynamic_obs_x_set[count_dynamic] = obs_x_set_gcs[i]
                dynamic_obs_y_set[count_dynamic] = obs_y_set_gcs[i]
                count_dynamic += 1

    #动态规划主函数
    def main_fun(self, obs_s_set, obs_l_set, plan_start_s, plan_start_l, plan_start_dl, plan_start_dll,
            w_cost_collision, w_cost_smooth_dl, w_cost_smooth_ddl, w_cost_smooth_dddl, w_cost_ref, row, col, sample_s, sample_l):
    #输入为： obs s l set 筛选后的障碍物s、l信息 ； plan_start s l dl ddl 规划起点信息 ； w_cost_obs, w_cost_smooth_dl, w_cost_smooth_ddl, w_cost_smooth_dddl, w_cost_ref 动态规划代价权重
    #           row, col 动态规划采样点的行数和列数（行数必须是奇数，列数必须为6） ； sample_s, sample_l 采样的s、l长度
    #输出为： dp_path_l, dp_path_S 动态规划的路径s、l，此路径不包含规划起点

        node_cost = np.ones((row, col)) * 9999      #代价矩阵，node_cost[i][j]代表从起点到第i行第j列的最小代价
        pre_node_index = np.zeros((row, col))        #路径矩阵，pre_node_index[i][j]代表从起点到i行j列的最优路径中前一个节点的行号
        
        #计算代价
        for i in range(row):
            node_cost[i][0] = self.CalcStartCost(plan_start_s, plan_start_l, plan_start_dl, plan_start_dll,i, sample_s, sample_l, w_cost_collision,
                                w_cost_smooth_dl, w_cost_smooth_ddl, w_cost_smooth_dddl, w_cost_ref, obs_s_set, obs_l_set, row)
        
        #main function
        for j in range(1, col):
            for i in range(row):
                cur_node_s = plan_start_s +  j * sample_s    #计算当前node[i][j]的s， l
                cur_node_l = ((row + 1)/2  - i) * sample_l
                for k in range(row):
                    pre_node_s =  plan_start_s + (j - 1) *sample_s      #计算前一列节点的s， l
                    pre_node_l = ((row + 1)/2  - k) * sample_l
                    cost_neighbour = self.CalcNeighbourCost(pre_node_s, pre_node_l, cur_node_s, cur_node_l, w_cost_collision,
                                    w_cost_smooth_dl, w_cost_smooth_ddl, w_cost_smooth_dddl, w_cost_ref, obs_s_set, obs_l_set)
                    
                    pre_min_cost = node_cost[k][j-1]    #起点到上一节点的最小代价为node_cost[k][j-1]
                    cost_temp = pre_min_cost + cost_neighbour    #起点到node_cost[i][j]的最小代价为node_cost[k][j-1] 再加上node[k][j-1]到node[i][j]的代价中最小的
                    if cost_temp < node_cost[i][j]:
                        node_cost[i][j] = cost_temp
                        pre_node_index[i][j] = k #记录上一列节点的行号

        index = 0   #找到node_cost最后一列中，cost最小的，记代价最小的行号为index
        min_cost = 9999
        for i in range(row):
            if node_cost[i][-1] < min_cost:
                min_cost = node_cost[i][-1]
                index = i
            #动态规划最优路径初始化
        dp_node_list_row  = np.zeros((col, 1))
        cur_index = index       #从后往前逆推
        for i in range(col):
            pre_index = pre_node_index[cur_index][len(pre_node_index[0]) - i]   #找到cur_Index前面节点的行号
            dp_node_list_row[col - i] = cur_index           #把cur_index放到dp_node_list_row（用于记录行号）对应的位置
            cur_index = pre_index               #赋值进行下一步
        
        dp_path_s = np.ones((15, 1)) * -1
        dp_path_l = np.ones((15, 1)) * -1
        for i in range(col):
            dp_path_s[i] = plan_start_s + i * sample_s
            dp_path_l[i] = ((row + 1) / 2 - dp_node_list_row[i]) * sample_l

        return dp_path_s, dp_path_l

        #计算相邻两个节点之间的cost
    def CalcNeighbourCost(self, pre_node_s, pre_node_l, cur_node_s, cur_node_l, w_cost_collision,
                                    w_cost_smooth_dl, w_cost_smooth_ddl, w_cost_smooth_dddl, w_cost_ref, obs_s_set, obs_l_set):
        start_l = pre_node_l
        start_dl,start_ddl = 0, 0
        end_l = cur_node_l
        end_dl, end_ddl = 0, 0
        start_s = pre_node_s
        end_s = cur_node_s
        coeff = self.CalcQuinticCoeffient(start_l, start_dl, start_ddl, end_l, end_dl, end_ddl, start_s, end_s)
        a0, a1, a2, a3, a4, a5 = coeff[0], coeff[1], coeff[2], coeff[3], coeff[4], coeff[5]
        ds, l, dl, ddl, dddl = np.zeros((10, 1)), np.zeros((10, 1)), np.zeros((10, 1)), np.zeros((10, 1)), np.zeros((10, 1))
        for i in range(len(ds)):
            ds[i] = start_s + (i - 1) * (end_s - start_s)/10
        l = a0 * np.ones((10, 1)) + a1 * ds + a2 * ds**2 + a3 * ds**3 + a4 * ds**4 + a5 * ds**5
        dl = a1 * np.ones((10, 1)) + 2 * a2 *ds + 3 * a3 * ds**2 + 4 * a4 * ds**3 + 5 * a5 * ds**4
        ddl = 2 * a2 * np.ones((10, 1)) + 6 * a4 * ds + 12 * a4 * ds **2 + 20 * a5 * ds**3
        dddl = 6 * np.ones((10, 1)) * a3 + 24 * a4*ds + 60 * a5 * ds**2
        cost_smooth = w_cost_smooth_dl * np.dot(dl.T, dl) + w_cost_smooth_ddl * np.dot(ddl.T, ddl) + w_cost_smooth_dddl * np.dot(dddl.T, dddl)
        cost_ref = w_cost_ref * np.dot(l.T, l)
        cost_collision = 0
        for i in range(len(obs_s_set)):
            if isnan(obs_s_set[i]):
                break
            #这里做了简化，认为障碍物为一个点
            dlon = np.ones((10, 1)) * obs_s_set[i] - ds
            dlat = np.ones((10, 1)) * obs_l_set[i] - l
            square_d = dlon**2 + dlat**2
            cost_collision_once = self.CalcObsCost(w_cost_collision, square_d)
            cost_collision = cost_collision + cost_collision_once
        cost = cost_collision + cost_smooth + cost_ref
        return cost


    #计算起点到第一层的cost
    def CalcStartCost(self, begin_s, begin_l, begin_dl, begin_ddl, cur_node_row, sample_s, sample_l, w_cost_collision,
                        w_cost_smooth_dl, w_cost_smooth_ddl, w_cost_smooth_dddl, w_cost_ref, obs_s_set, obs_l_set, row):
        start_l = begin_l
        start_dl = begin_dl
        start_ddl = begin_ddl
        end_l = ((row + 1) / 2 - cur_node_row) * sample_l
        end_dl, end_ddl = 0, 0
        start_s = begin_s
        end_s = begin_s + sample_s
        coeff = self.CalcQuinticCoeffient(start_l, start_dl, start_ddl, end_l, end_dl, end_ddl, start_s, end_s)
        a0, a1, a2, a3, a4, a5 = coeff[0], coeff[1], coeff[2], coeff[3], coeff[4], coeff[5]
        ds, l, dl, ddl, dddl = np.zeros((10, 1)), np.zeros((10, 1)), np.zeros((10, 1)), np.zeros((10, 1)), np.zeros((10, 1))
        for i in range(len(ds)):
            ds[i] = start_s + (i - 1) * sample_s/10
        l = a0 * np.ones((10, 1)) + a1 * ds + a2 * ds**2 + a3 * ds**3 + a4 * ds**4 + a5 * ds**5
        dl = a1 * np.ones((10, 1)) + 2 * a2 *ds + 3 * a3 * ds**2 + 4 * a4 * ds**3 + 5 * a5 * ds**4
        ddl = 2 * a2 * np.ones((10, 1)) + 6 * a4 * ds + 12 * a4 * ds **2 + 20 * a5 * ds**3
        dddl = 6 * np.ones((10, 1)) * a3 + 24 * a4*ds + 60 * a5 * ds**2
        cost_smooth = w_cost_smooth_dl * np.dot(dl.T, dl) + w_cost_smooth_ddl * np.dot(ddl.T, ddl) + w_cost_smooth_dddl * np.dot(dddl.T, dddl)
        cost_ref = w_cost_ref * np.dot(l.T, l)
        cost_collision = 0
        for i in range(len(obs_s_set)):
            if isnan(obs_s_set[i]):
                break
            #这里做了简化，认为障碍物为一个点
            dlon = np.ones((10, 1)) * obs_s_set[i] - ds
            dlat = np.ones((10, 1)) * obs_l_set[i] - l
            square_d = dlon**2 + dlat**2
            cost_collision_once = self.CalcObsCost(w_cost_collision, square_d)
            cost_collision = cost_collision + cost_collision_once

        cost = cost_collision + cost_smooth + cost_ref
        return cost
    
    #计算障碍物的碰撞代价,暂定超过4米的cost为零，4到3米的代价为1000/squard_d，在小于3米的cost为w_cost_collision
    def CalcObsCost(self, w_cost_collision, square_d):
        cost = 0 
        for i in range(len(square_d)):
            if square_d[i] < 16 and square_d[i] > 9:
                cost = cost + 1000/square_d[i]
            elif square_d < 9:
                cost = cost + w_cost_collision
        obs_cost = cost
        return obs_cost


    def CalcQuinticCoeffient(self, start_l, start_dl, start_ddl, end_l, end_dl, end_ddl, start_s, end_s):
        start_s2 = start_s * start_s
        start_s3 = start_s2 * start_s
        start_s4 = start_s3 * start_s
        start_s5 = start_s4 * start_s
        end_s2 = end_s * end_s
        end_s3 = end_s2 * end_s
        end_s4 = end_s3 * end_s
        end_s5 = end_s4 * end_s
        A = [[1, start_s,    start_s2,  start_s3,    start_s4,   start_s5],
            [0,     1,       2*start_s, 3*start_s2, 4*start_s3, 5*start_s4],
            [0,     0,      2,          6*start_s,  12*start_s2, 20*start_s3], 
            [1,     end_s,  end_s2,     end_s3,     end_s4,     end_s5],
            [0,     1,      2*end_s,    3*end_s2,   4*end_s3,   5*end_s4],
            [0,     0,      2,          6*end_s,    12*end_s2,  20*end_s3]]
        B = [[start_l], [start_dl], [start_ddl], [end_l], [end_dl], [end_ddl]]
        coeff = np.dot(np.matrix(A).I, np.matrix(B))
        return coeff
        
        #用于增密path，动态规划采样的行列未知，也就不知道动态规划计算的曲线，需要做缓冲
    def add_density(self, dp_path_s_init, dp_path_l_init, plan_start_s, plan_start_l, plan_start_dl, plan_start_ddl):
        ds = 1
        dp_path_s_final, dp_path_l_final, dp_path_dl_final, dp_path_ddl_final = np.zeros((60, 1)), np.zeros((60, 1)), np.zeros((60, 1)), np.zeros((60, 1))
        start_s, start_l = plan_start_s, plan_start_l   #设置五次多项式的起点
        start_dl, start_ddl = plan_start_dl, plan_start_ddl
        s_cur, l_temp, dl_temp, ddl_temp = [], [], [], []
        for i in range(len(dp_path_s_init)):
            if dp_path_s_init[i] == -1:
                break
            for j in range(10000):  #采样s
                s_node = start_s + ds * (j - 1)
                if s_node < dp_path_s_init[i]:
                    s_cur = [s_cur, s_node]
                else:
                    break

            end_s = dp_path_s_init[i]
            end_l = dp_path_l_init[i]
            end_dl, end_ddl = 0, 0
            coeff = self.CalcQuinticCoeffient(start_l, start_dl, start_ddl, end_l, end_dl, end_ddl, start_s, end_s)
            a0, a1, a2, a3, a4, a5 = coeff[0], coeff[1], coeff[2], coeff[3], coeff[4], coeff[5]
            l = a0 * np.ones((1, len(s_cur))) + a1 * s_cur + a2 * s_cur**2 + a3 * s_cur**3 + a4 * s_cur**4 + a5 * s_cur**5
            dl = a1 * np.ones((1, len(s_cur))) + 2 * a2 * s_cur + 3 * a3 * s_cur **2 + 4 * a4 * s_cur**3 + 5 * a5 *s_cur**4
            ddl = 2 * a2 * np.ones((1, len(s_cur))) + 6 * a3 * s_cur + 12 * a4 * s_cur**2 + 20 * a5 * s_cur**3
            l_temp = [l_temp, l]    #保存l， dl， ddl结果
            dl_temp = [dl_temp, dl]
            ddl_temp = [ddl_temp, ddl]

            satrt_s = end_s     #准备下一次循环
            start_l = end_l
            start_dl = end_dl
            start_ddl = end_ddl
            s_cur = []

        for k in range(len(l_temp)):
            if k == 60:
                break
            dp_path_s_final[k] = plan_start_s + ds * (k-1)
            dp_path_l_final[k] = l_temp[k]
            dp_path_dl_final[k] = dl_temp[k]
            dp_path_ddl_final[k] = ddl_temp[k]

        return dp_path_s_final, dp_path_l_final, dp_path_dl_final, dp_path_ddl_final
        
    #自然坐标系到直角坐标系
    def fernet2cartesian(self, s_set, l_set, dl_set, ddl_set, frenet_path_x, frenet_path_y, frenet_path_heading, frenet_path_kappa, index2s):
        x_set, y_set, heading_set, kappa_set = np.ones((128, 1)), np.ones((128, 1)), np.ones((128, 1)), np.ones((128, 1))
        for i in range(len(s_set)):
            if isnan(s_set[i]):
                break
            proj_x, proj_y, proj_heading, proj_kappa = self.CalcProjPoint(s_set[i], frenet_path_x, frenet_path_y, frenet_path_heading, frenet_path_kappa, index2s)
            nor = [-math.sin(proj_heading), math.cos(proj_heading)]
            point = [proj_x, proj_y] + l_set[i] * nor
            x_set[i] = point[0]
            y_set[i] = point[1]
            heading_set[i] = proj_heading + math.atan(dl_set[i] / (1 - proj_kappa * l_set[i]))
            kappa_set[i] = ((ddl_set[i] + proj_heading * dl_set[i] * math.tan(heading_set[i] - proj_heading)) * (math.cos(heading_set[i] - proj_heading))**2 / 
                                    (1 - proj_kappa * l_set[i]) + proj_kappa) * math.cos(heading_set[i] - proj_heading) / (1 - proj_kappa * l_set[i])

        return x_set, y_set, heading_set, kappa_set
    
    #计算在frenet坐标系下，点s， l在frenet坐标轴的投影的直角坐标（proj_x, proj_y, proj_heading, proj_kappa）
    def CalcProjPoint(self, s, frenet_path_x, frenet_path_y, frenet_path_heading, frenet_path_kappa, index2s):
        match_index = 0
        while self.index2s(match_index) < s:
            match_index += 1
        match_point = [frenet_path_x[match_index], frenet_path_y[match_index]]
        match_point_heading = frenet_path_heading[match_index]
        match_point_kappa = frenet_path_kappa[match_index]
        ds = s - self.index2s(match_index)
        match_tor = [math.cos(match_point_heading), math.sin(match_point_heading)]
        proj_point = match_point + ds * match_tor
        proj_heading = match_point_heading + ds * match_point_kappa
        proj_kappa = match_point_kappa
        proj_x, proj_y = proj_point[0], proj_point[1]
        return proj_x, proj_y, proj_heading, proj_kappa


    ###############
    #以下为速度规划部分(用于完善动态规划，不代表最终速度规划)
    ###############
    #计算trajectory的s 与x， y的对应关系，可以看作是trajectory index2s
    def tra_correspond_xy(self, trajectory_x, trajectory_y):
        n = len(trajectory_x)
        path_s = np.zeros((n, 1))
        sum = 0
        
        for i in range(1, len(trajectory_x)):
            if isnan(trajectory_x[i]):
                break
            sum += sqrt((trajectory_x[i] - trajectory_x[i-1])**2 + (trajectory_y[i] - trajectory_y[i-1])**2)
            path_s[i] = sum
        
        if i == n:
            path_s_end = path_s[-1]
        else:
            path_s_end = path_s[i - 1]

        return path_s_end, path_s

    #计算速度规划的初始条件
    def Calc_velocity_init(self, plan_start_vx, plan_start_vy, plan_start_ax, plan_start_ay, plan_start_heading):
        tor = [math.cos(plan_start_heading), math.sin(plan_start_heading)]
        v_t = np.dot(self.transposition(tor), [plan_start_vx, plan_start_vy])
        a_t = np.dot(self.transposition(tor), [plan_start_ax, plan_start_ay])
        plan_start_s_dot = v_t
        plan_start_s_dot2 = a_t
        return plan_start_s_dot, plan_start_s_dot2
    
    #简单的速度二次规划
    def velocity_qp(self, plan_start_s_dot, plan_start_s_dot2, s_end, recommend_T):
        # coder.extrinsic("quadprog")
        n = 51
        s_init, s_dot_init, s_dot2_init, relative_time_init = np.zeros((n, 1)), np.zeros((n, 1)), np.zeros((n, 1)), np.zeros((n, 1))
        Aeq, beq = np.zeros((3*n, 2*n -2)), np.zeros((2*n - 2, 1))
        lb = np.ones((3*n, 1))
        ub = lb
        dt = recommend_T/(n-1)
        A_sub = [[1, 0], [dt, 1], [(1/3)*dt**2, (0.5*dt)], [-1, 0], [0, -1], [(1/6)*dt**2, dt/2]]
        for i in range(n-1):
            for o in range(len(A_sub)):
                Aeq[3*i+o-2][2*i-1 : 2*i] = A_sub[o]
        
        for i in range(n):
            lb[3*i - 2] = -99999
            lb[3*i - 1] = 0
            lb[3*i] = -8
            ub[3*i - 2] = 99999
            ub[3*i - 1] = 50
            ub[3*i] = 2
        lb[0], lb[1], lb[2] = 0, plan_start_s_dot, plan_start_s_dot2
        ub[0], ub[1], ub[2] = lb[0], lb[1], lb[2]
        lb[3*n - 2] = s_end
        ub[3*n - 2] = s_end
        a3, a4 = np.zeros((3*n, 3*n))
        A4_sub = [0, 0, 1, 0, 0, -1]
        for i in range(n-1):
            a3[3*n][3*n] = 1
            for p in range(len(A4_sub)):
                a4[3*i+p-2][i:i] = A4_sub[p]
        a3[3*n][3*n] = 1
        hh = a3 + 100 * np.dot(a4, self.transposition(a4))
        hh = 2 * hh
        f = np.zeros((3 * n, 1))
        xx = quadprog(hh, f, [], [], self.transposition(Aeq), beq, lb, ub)
        for i in range(n):
            s_init[i] = xx[3*i - 2]
            s_dot_init[i] = xx[3*i - 1]
            s_dot2_init[i] = xx[3*i]
            relative_time_init[i] = (i -1)* dt
        
        return s_init, s_dot_init, s_dot2_init, relative_time_init

    #增密s， s_dot， s_dot2
    def add_density1(self, s_init, s_dot_init, s_dot2_init, relative_time_init):
        tt = relative_time_init[-1]
        n = 401
        dt = tt/(n-1)
        s, s_dot, s_dot2 = np.zeros((n, 1)), np.zeros((n, 1)), np.zeros((n, 1))
        for i in range(n):
            current_t = (i - 1)*dt
            for j in range(len(relative_time_init - 1)):
                if relative_time_init[j] <= current_t and relative_time_init[j+1] > current_t:
                    break
            x = current_t  - relative_time_init[j]
            s[i] = s_init[j] + s_dot_init[j] * x + (1/3) * s_dot2_init[j] * x**2 + (1/6) * s_dot2_init[j+1] * x**2
            s_dot[i] = s_dot_init[j] + 0.5*s_dot2_init[j]*x + 0.5 * s_dot2_init[j+1]*x
            s_dot2[i] = s_dot2_init[j] + (s_dot2_init[j+1] - s_dot2_init[j]) * x / (relative_time_init[j+2] - relative_time_init[j])
            # ................... 
        
        return s, s_dot, s_dot2, relative_time

    #合并path和speed
    def merge_path_velocity(self, s, s_dot, s_dot2, relative_time, current_time, path_s, trajectory_x_init, trajectory_y_init, trajectory_heading_init, trajectory):
        #path是60个点，speed有401个点，因此要做插值
        n = 401
        trajectory_x, trajectory_y, trajectory_heading, trajectory_kappa, trajectory_speed, trajectory_accel = np.zeros((n, 1)), np.zeros((n, 1)), np.zeros((n, 1)), np.zeros((n, 1)), np.zeros((n, 1)), np.zeros((n, 1))
        index = 1

        return trajectory_x, trajectory_y, trajectory_heading, trajectory_kappa, trajectory_speed, trajectory_accel, trajectory_time

    #轨迹拼接，将stitch trajectoruy拼到规划完成的轨迹上
    def stitch_trajectory(self, trajectory_x, trajectory_y, trajectory_heading, trajectory_kappa, trajectory_speed, trajectory_accel, trajectory_time, 
                        stitch_x, stitch_y, stitch_heading, stitch_kappa, stitch_speed, stitch_accel, stitch_time):
        n = 401 + 20
        trajectory_x_final, trajectory_y_final, trajectory_heading_final, trajectory_kappa_final, trajectory_speed_final, trajectory_accel_final, trajectory_time_final = np.zeros((n, 1)), np.zeros((n, 1)), np.zeros((n, 1)), np.zeros((n, 1)), np.zeros((n, 1)), np.zeros((n, 1)), np.zeros((n, 1))
        trajectory_x_final = [stitch_x, trajectory_x]
        trajectory_y_final = [stitch_y, trajectory_y]
        trajectory_heading_final = [stitch_heading, trajectory_heading]
        trajectory_kappa_final = [stitch_kappa, trajectory_kappa]
        trajectory_speed_final = [stitch_speed, trajectory_speed]
        trajectory_accel_final = [stitch_accel, trajectory_accel]
        trajectory_time_final = [stitch_time, trajectory_time]

        return trajectory_x_final, trajectory_y_final, trajectory_heading_final, trajectory_kappa_final, trajectory_speed_final, trajectory_accel_final, trajectory_time_final

    #######
    #control
    #######
    def trans2xryr(self, trajectory_x, trajectory_y, trajectory_heading, trajectory_kappa, trajectory_speed, trajectory_accel, trajectory_time, current_time):
        #该函数根据规划的轨迹计算出期望跟踪的点，由于控制一样有延迟， 所以控制发出的期望跟踪的点是current_time + 0.01
        control_time = current_time + 0.01

        #规划发出的轨迹有一部分是拼接的，在刚启动的时候，stitch_time = -1， 因此要把它去掉
        start_index = 1
        while trajectory_time[start_index] == -1:
            start_index += 1
        
        #首次运行时，规划执行完的时候控制已经执行十次了，那么这十次中是空轨迹，这样轨迹中的trajectory_time也全是零，会导致interp1插值报错
        if control_time > trajectory_time[start_index] and trajectory_time[start_index] != 0:
            xr = interp1(trajectory_time[start_index:], trajectory_x[start_index:], control_time)
            yr = interp1(trajectory_time[start_index:], trajectory_y[start_index:], control_time)
            thetar = interp1(trajectory_time[start_index:], trajectory_heading[start_index:], control_time)
            kappar = interp1(trajectory_time[start_index:], trajectory_kappa[start_index:], control_time)
            vr = interp1(trajectory_time[start_index:], trajectory_speed[start_index:], control_time)
            ar = interp1(trajectory_time[start_index:], trajectory_accel[start_index:], control_time)
        else:
            xr = trajectory_x[start_index]
            yr = trajectory_y[start_index]
            thetar = trajectory_heading[start_index]
            kappar = trajectory_kappa[start_index]
            vr = trajectory_speed[start_index]
            ar = trajectory_accel[start_index]

        return xr, yr, thetar, kappar, vr, ar