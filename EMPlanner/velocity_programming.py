import numpy as np
import math
from cmath import atan, isnan, nan, sqrt
import qpsolvers
from scipy import interpolate

from EMPlanner.quadratic_programming import qp

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
    def speed_dp(self, obs_st_s_in_set, obs_st_s_out_set, obs_st_t_in_set, obs_st_t_out_set, reference_speed, w_cost_ref_speed, w_cost_accel, w_cost_obs, plan_start_s_dot):
        #时间从0到8开始规划，最多8秒 ； s的范围从0开始到路径规划的path的总长度 ； 为减少算力，采用非均匀采样，越小的越密，越大的越稀疏
        s_list = [np.arange(start=0, stop=4.5, step=0.5), np.arange(start=5.5, stop=14.5, step=1), np.arange(start=16, stop=29.5, step=1.5), np.arange(start=32, stop=54.5, step=2.5)]
        s_list = np.reshape(s_list, (len(s_list, s_list[0])))
        t_list = np.arange(start=0.5, stop=8, step=0.5)
            #s采样40个点，t采样16个点
        dp_st_cost = np.ones((40, 16)) * 99999  #代价矩阵
        dp_st_node = np.ones((40, 16))
        dp_st_s_dot = np.zeros((40, 16))        #表示从起点到i，j点的最优路径的末速度

        for i in range(len(s_list)):    #计算从dp起点到第一列的cost
            dp_st_cost[i][0] = self.CalcDPCost(0, 0, i, 1, obs_st_s_in_set, obs_st_s_out_set, obs_st_t_in_set, obs_st_t_out_set, w_cost_ref_speed, reference_speed, w_cost_accel, w_cost_obs, plan_start_s_dot, s_list, t_list, dp_st_s_dot)
            s_end, t_end = self.CalcSTCoordinate(i, 1, s_list, t_list)
            dp_st_s_dot[i][0] = s_end / t_end
        
        for i in range(1, len(t_list)):
            for j in range(len(s_list)):
                cur_row = j
                cur_col = i
                for k in range(len(s_list)):
                    pre_row = k
                    pre_col = i - 1
                    cost_temp = self.CalcDPCost(self, pre_row, pre_col, cur_row, cur_col, obs_st_s_in_set, obs_st_s_out_set, obs_st_t_in_set, obs_st_t_out_set, w_cost_ref_speed, reference_speed, w_cost_accel, w_cost_obs, plan_start_s_dot, s_list, t_list, dp_st_s_dot)

                    if cost_temp + dp_st_cost[pre_row][pre_col] < dp_st_cost[cur_row][cur_col]:
                        dp_st_cost[cur_row][cur_col] = cost_temp + dp_st_cost[pre_row, pre_col]
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
                dis3 = abs(vector1[0] * vector3[1] - vector1[1] * vector3[0]) / sqrt(np.dot(vector3.T, vector3))
                if (np.dot(vector1.T, vector3) > 0 and np.dot(vector2.T, vector3) > 0) or (np.dot(vector1.T, vector3) < 0 and np.dot(vector2.T,vector3) < 0):
                    min_dis = min(dis1, dis2)
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

        #计算出s， s_dot, 的上下界，开辟凸空间供二次规划使用
    def Convex_space(self, dp_speed_s, dp_speed_t, path_index2s, obs_st_s_in_set, obs_st_s_out_set, obs_st_t_in_set, obs_st_t_out_set, trajectory_kappa_init, max_lateral_accel):
        n = 16  #dp_speed_s t 最多有16个
        s_lb, s_ub = np.ones((n, 1)) * -99999, np.ones((n, 1)) * 999999
        s_dot_lb, s_dot_ub = np.ones((n, 1)) * -99999, np.ones((n, 1)) * 99999
        path_index2s_end_index = len(path_index2s)
        dp_speed_end_index = len(dp_speed_s)

        for k in range(1, len(path_index2s)):   #path_index2s也是有缓冲的，找到有效的path_index2s的末尾的位置，赋值给path_index2s_end_index
            if path_index2s[k] == 0 and path_index2s[k-1] !=0:
                path_index2s_end_index = k - 1
                break
            path_index2s_end_index = k

        for k in range(len(dp_speed_s)):    #找到dp_speed_s中有效的位置，赋值给dp_speed_end_index
            if isnan(dp_speed_s):
                dp_speed_end_index = k - 1
                break

        for i in range(n):  #施加车辆动力学约束
            if isnan(dp_speed_s[i]):
                break
            cur_s = dp_speed_s[i]
            cur_kappa = interpolate(path_index2s[0:path_index2s_end_index], trajectory_kappa_init[0:path_index2s_end_index], cur_s) #插值找cur_s所对应的曲率
            max_speed = sqrt(max_lateral_accel / abs(cur_kappa) + 1e-10)    #由a = v**2 * k
            min_speed = 0
            s_dot_lb[i] = min_speed
            s_dot_ub[i] = max_speed

        for i in range(len(obs_st_s_in_set)):
            if isnan(obs_st_s_in_set[i]):
                continue
            obs_t = (obs_st_t_in_set[i] + obs_st_t_out_set[i]) / 2  #取s t直线的中点， 作为obs_s obs_t的坐标
            obs_s = (obs_st_s_in_set[i] + obs_st_s_out_set[i]) / 2
            obs_speed = (obs_st_s_out_set[i] - obs_st_s_in_set[i]) / (obs_st_t_out_set[i] - obs_st_t_in_set[i])     #计算障碍物的纵向速度
            dp_s = interpolate([0, dp_speed_t[1:dp_speed_end_index]],  [0, dp_speed_s[1:dp_speed_end_index], obs_t])    #插值找到当t = obs_t时， dp_speed_s 的值
            for j in range(len(dp_speed_t) - 1):
                if dp_speed_t[0] > obs_st_t_in_set[i]:  #如果障碍物切入时间比0.5秒还要短，那么t_lb_index = 1
                    break
                elif dp_speed_t[j] <= obs_st_t_in_set[i] and dp_speed_t[j+1]:   #否则遍历dp_speed_t 找到与obs_st_t_in_set[i]最近的点的编号
                    break

            t_lb_index = j   #将点的编号赋值给j
            for j in range(len(dp_speed_t) - 1):    #找到dp_speed_t中与obs_st_t_out_set[i]最近的时间，并将此时间的编号赋值给t_ub_index
                if dp_speed_t[0] > obs_st_t_out_set[i]:
                    break
                elif dp_speed_t[j] <= obs_st_t_out_set[i] and dp_speed_t[j + 1] > obs_st_t_out_set[i]:
                    break
            
            t_ub_index = j  #这里稍微做个缓冲， 把t_lb_index稍微缩小一些，t_ub_index稍微放大一些
            t_lb_index = max(t_lb_index -2 , 3) #最低为3， 因为碰瓷没法处理
            t_ub_index = min(t_ub_index + 2, dp_speed_end_index)
            if obs_s > dp_s:    #决策为减速避让
                for m in range(t_lb_index, t_ub_index+1):     #在t_lb_index, t_ub_index的区间上，s的上界不可以超过障碍物st直线
                    dp_t = dp_speed_t[m]
                    s_ub[m] = min(s_ub[m], obs_st_s_in_set[i] + obs_speed * (dp_t - obs_st_t_in_set[i]))

            else:   #决策为加速超车
                for m in range(t_lb_index, t_ub_index+1):
                    dp_t = dp_speed_t[m]
                    s_lb[m] = max(s_lb[m], obs_st_s_in_set[i] + obs_speed * (dp_t - obs_st_t_in_set[i])) 
        
        return s_lb, s_ub, s_dot_lb, s_dot_ub

    #速度二次规划
    def qp_velocity(self, plan_start_s_dot, plan_start_s_dot2, dp_speed_s, dp_speed_t, s_lb, s_ub, s_dot_lb, s_dot_ub, w_cost_s_dot2, w_cost_v_ref, w_cost_jerk, speed_reference):
        #输入：规划起点plan_start_s_dot, plan_start_s_dot2 ； 动态规划结果：dp_speed_s, dp_speed_t ； 凸空间约束s_lb, s_ub, s_dot_lb, s_dot_ub ； 加速度、推荐速度、jerk代价权重w_cost_s_dot2, w_cost_v_ref, w_cost_jerk
        #输出：速度曲线
        # coder.extrinsic("quadprog")
        dp_speed_end = 16   #由于dp的结果未必是16， 该算法将计算dp_speed_end到底是多少
        for i in range(len(dp_speed_s)):
            if isnan(dp_speed_s[i]):
                dp_speed_end = i - 1
                break
        
        n = 17 #由于dp_speed_end实际上是不确定的，但是输出要求必须是确定长度的值，因此输出初始化选择dp_speed_end的最大值 + 规划起点作为输出初始化的规模
        qp_s_init, qp_s_dot_init, qp_s_dot2_init, relative_time_init = np.ones((n, 1)) * nan, np.ones((n, 1)) * nan, np.ones((n, 1)) * nan, np.ones((n, 1)) * nan
        s_end = dp_speed_s[dp_speed_end]
        recommend_T = dp_speed_t[dp_speed_end]  #此时dp_speed_end表示真正有效的dp_speed_t的元素个数，取出dp_speed_t有效的最后一个元素作为规划的时间终点，记为recommend_T
        qp_size = dp_speed_end + 1      #qp的规模应该是dp的有效元素的个数 + 规划起点
        Aeq, beq = np.zeros((3*qp_size, 2*qp_size - 2)), np.zeros((2*qp_size - 2, 1))
        lb, ub = np.ones((3*qp_size, 1)), lb
        dt = recommend_T / dp_speed_end
        A_sub = [[1, 0],
                [dt, 1],
                [(1/3)*dt**2, (1/2)*dt],
                [-1, 0],
                [0, -1], 
                [(1/6)*dt**2, dt/2]]
        
        for i in range(qp_size - 1):
           Aeq[3*i:3*i+5, 2*i:2*i+2] = A_sub
        
        A, b = np.zeros((qp_size-1, 3*qp_size)), np.zeros((qp_size - 1, 1))     #不允许倒车约束，即s(i) - s(i+1) <= 0
        for i in range(qp_size - 1):
            A[i, 3*i], A[i, 3*i-1] = 1, -1
        
        for i in range((1, qp_size)):   #由于生成的凸空间约束s_lb s_ub 不带起点(动态规划不带起点而二次规划带起点），所以lb(i) = s_lb(i-1)， 以此类推, 基于车辆动力学，最小加速度为-6，最大加速度为4
            lb[3*i] = s_lb[i]
            lb[3*i + 1] = s_dot_lb[i]
            lb[3*i + 2] = -6
            ub[3*i] = s_ub[i]
            ub[3*i + 1] = s_dot_ub[i]
            ub[3*i + 2] = 4
        
        lb[0], lb[1], lb[2] = 0, plan_start_s_dot, plan_start_s_dot2    #起点约束
        ub[0], ub[1], ub[2] = lb[0], lb[1], lb[2]

        A_s_dot2, A_jerk, A_ref = np.zeros((3*qp_size, 3*qp_size)), np.zeros((3*qp_size, qp_size - 1)), np.zeros((3*qp_size, 3*qp_size))
        A4_sub = [0, 0, 1, 0, 0, -1]
        for i in range(qp_size):
            A_s_dot2[3*i + 2][3*i + 2], A_ref[3*i + 1][3*i + 1] = 1, 1
        for i in range(qp_size-1):
            A_jerk[3*i:3*i+5, i:i] = A4_sub
        H = w_cost_s_dot2 * np.dot(A_s_dot2, A_s_dot2.T) + w_cost_jerk * np.dot(A_jerk, A_jerk.T) + w_cost_v_ref * np.dot(A_ref, A_ref.T)
        H = 2 * H
        f = np.zeros((3*qp_size, 1))
        for i in range(qp_size):
            f[3*i+1] = -2 * w_cost_v_ref * speed_reference

        X = qpsolvers.solve_qp(H, f, A, b , Aeq, beq, lb, ub)

        for i in range(qp_size):
            qp_s_init[i] = X[3*i]
            qp_s_dot_init[i] = X[3*i + 1]
            qp_s_dot2_init[i] = X[3*i + 2]
            relative_time_init[i] = (i - 1) * dt
        
        return qp_s_init, qp_s_dot_init, qp_s_dot2_init, relative_time_init

    # 该函数将增密s s_dot s_dot2
    def add_density_s(self, s_init,s_dot_init,s_dot2_init,relative_time_init):
        t_end  = len(relative_time_init)
        for i in range(t_end):
            if isnan(relative_time_init[i]):
                t_end = i - 1
                break
        T = relative_time_init[t_end]
        n = 401
        dt = T / (n - 1)
        s = np.zeros(n, 1)
        s_dot = np.zeros(n, 1)
        s_dot2 = np.zeros(n, 1)
        relative_time = np.zeros(n, 1)

        for i in range(n):
            current_t = (i - 1) * dt
            for j in range(t_end - 1):
                if relative_time_init[j] <= current_t and relative_time_init[j+1] > current_t:
                    break
            
            x = current_t - relative_time_init[j]
            s[i] = s_init[i] + s_dot_init[j] * x + (1/3) * s_dot2_init[j] * x**2 + (1/6) * s_dot2_init[j+1] * x**2
            s_dot[i] = s_dot_init[j] + 0.5 * s_dot2_init[j] * x + 0.5 * s_dot2_init[j+1] * x
            s_dot2[i] = s_dot2_init[j] + (s_dot2_init[j+1] - s_dot2_init[j]) * x / (relative_time_init[j+1] - relative_time_init[j])
            relative_time[i] = current_t

        return s,s_dot,s_dot2,relative_time

        # 该函数将合并path和speed
        # 由于path 是 60个点，speed 有 401个点，合并后，path和speed有401个点，因此需要做插值
    def merge_path_velocity(self, s,s_dot,s_dot2,relative_time,current_time,path_s,trajectory_x_init,trajectory_y_init,trajectory_heading_init,trajectory_kappa_init):
        n = 401
        trajectory_x = np.zeros(n, 1)
        trajectory_y = np.zeros(n, 1)
        trajectory_heading  = np.zeros(n, 1)
        trajectory_kappa = np.zeros(n, 1)
        trajectory_speed = np.zeros(n, 1)
        trajectory_accel = np.zeros(n, 1)
        index = 1

        while not isnan(trajectory_x_init[index]):
            index += 1
        index -= 1
        trajectory_time = np.zeros(n, 1)

        for i in range(n-1):
            trajectory_x[i] = interpolate(path_s[1:index], trajectory_x_init[1:index], s[i])
            trajectory_y[i] = interpolate(path_s[1:index], trajectory_y_init[1:index], s[i])
            trajectory_heading[i] = interpolate(path_s[1:index], trajectory_heading_init[1:index], s[i])
            trajectory_kappa[i] = interpolate(path_s[1:index], trajectory_kappa_init[1:index], s[i])
            trajectory_time[i] = relative_time[i] + current_time
            trajectory_speed[i] = s_dot[i]
            trajectory_accel[i] = s_dot2[i]
        
        trajectory_x[-1] = trajectory_x_init[-1]
        trajectory_y[-1] = trajectory_y_init[-1]
        trajectory_heading[-1] = trajectory_heading_init[-1]
        trajectory_kappa[-1] = trajectory_kappa_init[-1]
        trajectory_time[-1] = relative_time[-1] + current_time
        trajectory_speed[-1] = s_dot[-1]
        trajectory_accel[-1] = s_dot2[-1]
        
        return trajectory_x,trajectory_y,trajectory_heading,trajectory_kappa,trajectory_speed,trajectory_accel,trajectory_time

    def stitch_trajectory(self, trajectory_x, trajectory_y,trajectory_heading,trajectory_kappa,trajectory_speed,trajectory_accel,trajectory_time, stitch_x,stitch_y,stitch_heading,stitch_kappa,stitch_speed,stitch_accel,stitch_time):
        n = 401 + 20
        trajectory_x_final = np.zeros(n, 1)
        trajectory_y_final = np.zeros(n, 1)
        trajectory_heading_final = np.zeros(n, 1)
        trajectory_kappa_final = np.zeros(n, 1)
        trajectory_speed_final = np.zeros(n, 1)
        trajectory_accel_final = np.zeros(n, 1)
        trajectory_time_final = np.zeros(n, 1)

        trajectory_x_final = [stitch_x, trajectory_x]
        trajectory_y_final = [stitch_y, trajectory_y]
        trajectory_heading_final = [stitch_heading, trajectory_heading]
        trajectory_kappa_final = [stitch_kappa, trajectory_kappa]
        trajectory_speed_final = [stitch_speed, trajectory_speed]
        trajectory_accel_final = [stitch_accel, trajectory_accel]
        trajectory_time_final = [stitch_time, trajectory_time]
        return trajectory_x_final,trajectory_y_final,trajectory_heading_final,trajectory_kappa_final, trajectory_speed_final,trajectory_accel_final,trajectory_time_final