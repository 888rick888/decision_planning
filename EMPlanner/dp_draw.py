import numpy as np
import matplotlib.pyplot as plt
import random



class dp_self:
    def __init__(self) -> None:
        self.sample_s = 1
        self.sample_l = 1
        self.row = 9
        self.col = 6

    def calccost(self):
        return random.random()

    def dp(self):
        #代价矩阵记录从起点到该点的最小代价， 点索引矩阵记录到该点的最优路径中的前一个点的行号
        node_cost, pre_node_index = np.ones((self.row, self.col)) * 99999, np.zeros((self.row, self.col))
        for i in range(self.row):
            node_cost[i][0] = self.calccost()
        
        for j in range(1, self.col):
            for i in range(self.row):
                cur_node_s = j * self.sample_s
                cur_node_l = ((self.row + 1) / 2 - i) * self.sample_l
                for k in range(self.row):
                    pre_node_s = (j - 1) * self.sample_s
                    pre_node_l = ((self.row + 1) / 2 - k) * self.sample_l
                    cost_neighbour = self.calccost()
                    pre_min_cost = node_cost[k][j-1]
                    cost_temp = pre_min_cost + cost_neighbour
                    if cost_temp < node_cost[i][j]:
                        node_cost[i][j] = cost_temp
                        pre_node_index[i][j] = k
        
        index = 0
        min_cost = 9999
        path_s, path_l = np.zeros((self.row, 1)), np.zeros((self.row, 1))
        dp_node_list_row = np.zeros((self.col, 1))  #记录该列的最优点的行号
        for p in range(self.row):
            if node_cost[p][-1] < min_cost:
                min_cost = node_cost[p][-1]
                index = p

        cur_index = index
        for i in range(self.col):
            pre_index  = pre_node_index[int(cur_index)][len(pre_node_index[0]) - i - 1]
            dp_node_list_row[self.col - i -1] = cur_index
            cur_index = pre_index
        for i in range(self.col):
            path_s[i] = i * self.sample_s
            path_l[i] = ((self.row + 1) / 2 - dp_node_list_row[i]) * self.sample_l
        
        return path_s, path_l

    def draw(self, point, path_s, path_l):
        plt.figure()
        plt.title("dynamic planning")

        point_x = [i for i in range(len(point))]
        point_y = [i for i in range(len(point))]
        # plt.scatter(point_x, point_y)
        # plt.plot(path_s, path_l)
        plt.plot(path_l)
        plt.show()

    def main(self):
        path_s, path_l = self.dp()
        point = np.zeros((self.row, self.col))
        self.draw(point, path_s, path_l)



if __name__=="__main__":
    dp = dp_self()
    dp.main() 