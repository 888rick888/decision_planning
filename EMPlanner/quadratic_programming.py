import numpy as np
import math
from cmath import atan, isnan, sqrt
import cvxopt
import qpsolvers

W_COST_L = 5
W_COST_DL = 50000
W_COST_DDL = 300
W_COST_DDDL = 15
W_COST_CENTRE = 1200
W_COST_END_L = 40
W_COST_END_DL = 40
W_COST_END_DDL = 40

class qp:
    def __init__(self) -> None:
        pass
    #该函数将输出每个离散的dp_path_s中的点s所对应的l的边界l_min,l_max
    def get_boundary(self, dp_path_s, dp_path_l, static_obs_s_set, static_obs_l_set, static_obs_length, static_obs_width):
        #输入动态规划的曲线dp_path_s, dp_path_l， 障碍物中心点的坐标static_obs_s_set, static_obs_l_set; 障碍物的长宽static_obs_length, static_obs_width
        #一般真实障碍物投影到frenet后，长宽会扭曲， 这里近似用直角坐标系的static_obs_length, static_obs_width的值代替frenet坐标下的值， 所以这套算法不能处理参考线曲率过大的场景（需要r>=200）， 大曲率的场景，需要对参考线模块、障碍物投影模块、速度规划模块都要做很多特殊处理

        l_min, l_max = np.ones((60, 1)) * -6, np.ones((60, 1)) * 6  #如果无障碍，l的边界一般默认为±6

        for i in range(len(static_obs_s_set)):
            if isnan(static_obs_s_set):
                break

            obs_s_min = static_obs_s_set[i] - static_obs_length / 2     #计算障碍物头尾部的s
            obs_s_max = static_obs_s_set[i] + static_obs_length / 2

            start_index = self.find_near_index(dp_path_s, obs_s_min)    #计算障碍物在dp_path_s上的位置
            end_index = self.find_near_index(dp_path_s, obs_s_max)

            centre_index = self.find_near_index(dp_path_s, static_obs_s_set[i])

            if start_index == 1 and end_index == 1:
                continue
            path_l = dp_path_l[centre_index]
            if path_l > static_obs_l_set[i]:    #l大于障碍物中心点，则向左绕过障碍物
                for j in range(start_index, end_index + 1):
                    l_min[j] = max(l_min[j], l_min[j] + static_obs_l_set[i] + static_obs_width / 2)
            else:
                for j in range(start_index, end_index + 1):
                    l_max[j] = min(l_max[j], static_obs_l_set[j] - static_obs_width / 2)


        return l_min, l_max

    #该函数将找到dp_path_s上所有点中， 与pns_S最近的点，并返回该点在dp_path_s的编号
    def find_near_index(self, dp_path_s, obs_s):
        #dp_path_s是动态规划的结果，值在0到59， 但是obs_s是障碍物在参考线上的投影点，所以obs_s有可能小于零，也有可能大于59
        if dp_path_s[0] >= obs_s:
            y = 0
        elif dp_path_s[-1] < obs_s:
            y = 60
        else:
            index = 0
            while dp_path_s[index] < obs_s:
                index += 1

        if dp_path_s[index] - obs_s > obs_s - dp_path_s[index - 1]:
            y = index - 1
        else:
            y = index

        return y

    #路径二次规划   0.5 * x' Hx + f' *x = min  , subject to A * x <= b , Aeq * x = beq , lb <= x <= ub
    def mian_fun(self, l_min, l_max, w_cost_l, w_cost_dl, w_cost_ddl, w_cost_dddl, w_cost_centre, w_cost_end_l, w_cost_end_dl, w_cost_end_ddl,
                host_d1, host_d2, host_w, plan_start_l, plan_start_dl, plan_start_ddl, dp_path_l_final):
        #输入为：l_min, l_max 点的凸空间 ； w_cost_l 参考线代价，w_cost_dl ddl dddl 光滑性代价， w_cost_centre凸空间中央代价
        # w_cost_end_l dl ddl 终点的状态代价， host_d1 d2 车辆质心到前后轴的距离， host_w 车的宽度， plan_start 规划起点
        #输出：qp_path_l dl ddl 二次规划输出曲线

        n = 60
        qp_path_l, qp_path_dl, qp_path_ddl = np.ones((n, 1)), np.ones((n, 1)), np.ones((n, 1))
        H_L, H_DL, H_DDL, H_DDDL, H_L_END, H_DL_END, H_DDL_END = np.zeros((3*n, 3*n)), np.zeros((3*n, 3*n)), np.zeros((3*n, 3*n)), np.zeros((n-1, 3*n)), np.zeros((3*n, 3*n)), np.zeros((3*n, 3*n)), np.zeros((3*n, 3*n))
        Aeq, beq, A, b = np.zeros((2*n-2, 3*n)), np.zeros((2*n-2, 1)), np.zeros((8*n, 3*n)), np.zeros((8*n, 1))
        end_l_desire, end_dl_desire, end_ddl_desire = 0, 0, 0

        ds = 1  #纵向间隔
        Aeq_sub = [[1, ds, ds**2/3, -1, 0, ds**2/6], 
                    [0, 1, ds/2, 0, -1, ds/2]]
        
        d1, d2, w = host_d1, host_d2, host_w
        A_sub = [[1, d1, 0],
                [1, d1, 0],
                [1, -d2, 0],
                [1, -d2, 0],
                [-1, -d1, 0],
                [-1, -d1, 0],
                [-1, d2, 0],
                [-1, d2, 0]]
        for i in range(n - 1):
            row = 2*i - 1 +2
            col = 3*i - 2 +3
            Aeq = np.matrix(Aeq)
            Aeq[row:row+2, col:col+6] = Aeq_sub

        for i in range(n):
        # for i in range(1， n):    #二次规划奔溃的解决方法之二
            row = 8*i - 7 +8
            col = 3*i - 2 +3
            A = np.matrix(A)
            A[row:row + 8, col:col + 3] = A_sub

        front_index = math.ceil(d1/ds)
        back_index = math.ceil(d2/ds)

        for i in range(n):
        # for i in range(1， n):    #二次规划奔溃的解决方法之二
            index1 = min(i + front_index, n)    #左前右前
            index2 = max(i - back_index, 1)     #左后右后
            b[8*i - 7:8*i+1 , 0] = [[l_max[index1] - w/2],[l_max[index1] + w/2], [l_max[index2] - w/2], [l_max[index2] + w/2], [-l_min[index1] + w/2], [-l_min[index1] - w/2], [-l_min[index2] + w/2], [-l_min[index2] - w/2]]

        lb = np.ones((3*n, 1)) * -99999
        ub = np.ones((3*n, 1)) * 99999
        lb[0], lb[1], lb[2] = plan_start_l, plan_start_dl, plan_start_ddl
        ub[0], ub[1], ub[2] = lb[0], lb[1], lb[2]

        for i in range(n):
            H_L[3*i - 2, 3*i - 2], H_DL[3*i - 1, 3*i - 1], H_DDL[3*i, 3*i] = 1, 1, 1
        
        H_CENTRE = H_L
        H_dddl_sub = [0, 0, 1, 0, 0, -1]
        for i in range(n-1):
            row = i
            col = 3*i - 2
            H_DDDL[row, col:col+6] = H_dddl_sub

        H_L_END[3*n - 2, 3*n - 2], H_DL_END[3*n - 1, 3*n -1], H_DDL_END[3*n, 3*n] = 1, 1, 1
        H = w_cost_l * np.dot(H_L.T, H_L) + w_cost_dl * np.dot(H_DL.T, H_DL) + w_cost_ddl * np.dot(H_DDL.T, H_DDL) + w_cost_dddl * np.dot(H_DDDL.T, H_DDDL) + w_cost_centre * np.dot(H_CENTRE.T, H_CENTRE) + w_cost_end_l * np.dot(H_L_END.T, H_L_END) + w_cost_end_dl * np.dot(H_DL_END.T, H_DL_END) + w_cost_ddl * np.dot(H_DDL_END.T, H_DDL_END)
        H = 2 * H

        f = np.zeros((3*n, 1))
        centre_line = 0.5 * (l_min + l_max)
        # centre_line = dp_path_l_final     #二次规划奔溃的解决方法之一

        for i in range(n):
            f[3*i -2] = -2 * centre_line[i]
        f =  w_cost_centre * f

        f[3*n - 2] = f[3*n - 2] - 2 * end_l_desire * w_cost_end_l 
        f[3*n - 1] = f[3*n - 1] - 2 * end_dl_desire * w_cost_end_dl
        f[3*n] = f[3*n] - 2 * end_ddl_desire * w_cost_end_ddl

        # X = cvxopt.solvers.qp(H, f, A, b , Aeq, beq, lb, ub)  #使用这个库需要将lb和ub放到A和b中
        X = qpsolvers.solve_qp(H, f, A, b , Aeq, beq, lb, ub)
        #缺乏曲率的约束，可以用|dl[i+1] - dl[i]|/ds <= kappa_max 近似约束曲率
        for i in range(n):
            qp_path_l[i] = X[3*i - 2]
            qp_path_dl[i] = X[3*i - 1]
            qp_path_ddl[i] = X[3*i]

        return qp_path_l, qp_path_dl, qp_path_ddl
