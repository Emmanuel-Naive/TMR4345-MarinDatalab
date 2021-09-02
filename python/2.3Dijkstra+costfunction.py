import matplotlib.pyplot as plt
import numpy as np
import math


class Dijkstra:

    def __init__(self, ox, oy,cosfun,sx,sy):
        # 初始化地图的情况
        self.min_x = None
        self.max_x = None
        self.min_y = None
        self.max_y = None
        self.x_grid_num = None
        self.y_grid_num = None
        self.obstacle_map = None

        self.calc_obstacle_grid_map(ox, oy)         # 构建环境栅格地图
        self.robot_motion = self.get_motion_model(sx,sy,cosfun)

    def calc_obstacle_grid_map(self, ox, oy):
        """ 构建环境栅格地图 """
        # 1. 获取环境的 上、 下、 左、 右 四个边界值
        self.min_x = round(min(ox))
        self.max_x = round(max(ox))
        self.min_y = round(min(oy))
        self.max_y = round(max(oy))

        # 2. 根据四个边界值和栅格的大小计算 x, y 方向上 栅格的数量
        self.x_grid_num = round(self.max_x - self.min_x)
        self.y_grid_num = round(self.max_y - self.min_y)

        # 3. 初始化环境栅格地图
        self.obstacle_map = [[False for _ in range(self.x_grid_num)] for _ in range(self.y_grid_num)]

    def planning(self, sx, sy, gx, gy,cosfun):
        """ 进行路径规划 """

        # 1. 将机器人的坐标进行结点化
        sx_index = self.calc_xy_index(sx, self.min_x)
        sy_index = self.calc_xy_index(sy, self.min_y)
        gx_index = self.calc_xy_index(gx, self.min_x)
        gy_index = self.calc_xy_index(gy, self.min_y)
        start_node = self.Node(sx_index, sy_index, 0.0, -1)
        goal_node = self.Node(gx_index, gy_index, 0.0, -1)

        # 2. 初始化 open_set, close_set,并将起点放进 open_set 中
        open_set, close_set = dict(), dict()
        open_set[self.calc_index(start_node)] = start_node

        # 3.开始循环
        while True:
            # (1). 取 open_set 中 cost 最小的结点
            c_id = min(open_set, key=lambda o: open_set[o].cost)
            current = open_set[c_id]
           
            if show:  # pragma: no cover
                plt.plot(self.calc_position(current.x, self.min_x),
                         self.calc_position(current.y, self.min_y), "xc")
                # for stopping simulation with the esc key.
                # plt.gcf().canvas.mpl_connect(
                #     'key_release_event',
                #     lambda event: [exit(0) if event.key == 'escape' else None])
                if len(close_set.keys()) % 10 == 0:
                    plt.pause(0.001)

            # (2). 判断该节点是否为终点
            if current.x == goal_node.x and current.y == goal_node.y:
                print('Find Goal!')
                goal_node.parent_index = current.parent_index
                goal_node.cost = current.cost
                break

            # (3). 将该节点从 open_set 中取出，并加入到 close_set 中
            del open_set[c_id]
            close_set[c_id] = current

            # (4). 根据机器人的运动模式，在栅格地图中探索当前位置出发到达的下一可能位置
            self.robot_motion = self.get_motion_model(current.x,current.y,cosfun)
            for move_x, move_y, move_cost in self.robot_motion:
                node = self.Node(current.x + move_x,
                                 current.y + move_y,
                                 current.cost + move_cost, c_id)
                n_id = self.calc_index(node)

                if n_id in close_set:
                    continue

                if not self.verify_node(node):
                    continue

                if n_id not in open_set:
                    open_set[n_id] = node   # 发现新的结点
                else:
                    if open_set[n_id].cost >= node.cost:
                        # 当前节点的路径到目前来说是最优的，进行更新
                        open_set[n_id] = node

        rx, ry = self.calc_final_path(goal_node, close_set)

        return rx, ry

    def calc_final_path(self, goal_node, close_set):
        """ 从终点开始进行回溯，生成从起点到终点的最优路径 """
        rx = [self.calc_position(goal_node.x, self.min_x)]
        ry = [self.calc_position(goal_node.y, self.min_y)]
        
        parent_index = goal_node.parent_index
        while parent_index != -1:
            n = close_set[parent_index]
            rx.append(self.calc_position(n.x, self.min_x))
            ry.append(self.calc_position(n.y, self.min_y))
            parent_index = n.parent_index
            
        return rx, ry

    class Node:
        def __init__(self, x, y, cost,parent_index):
            self.x = x      # 栅格的 x 轴索引
            self.y = y      # 栅格的 y 轴索引
            self.cost = cost       # g(n)
            self.parent_index = parent_index       # 当前节点的父节点
        #
        # def __str__(self):
        #     return str(self.x) + "," + str(self.y) + "," + str(self.cost) + "," + str(self.parent_index)

    def calc_index(self, node):
        """
        将栅格化后的地图进行编号索引，从左下角向右一行一行进行编号索引,如下面示例
        [7, 8, 9]
        [4, 5, 6]
        [1, 2, 3]
        """
        index = node.y * self.x_grid_num + node.x
        return index

    def calc_xy_index(self, pos, min_p):
        """ 将机器人在二维环境地图中的坐标转化成栅格地图中的坐标 """
        index = round(pos - min_p)
        return index

    def calc_position(self, index, min_p):
        """ 将栅格地图的坐标转化成在真实环境中的坐标 """
        pos = min_p + index
        return pos

    def verify_node(self, node):
        """ 验证机器人的当前位置是否合理 """
        px = self.calc_position(node.x, self.min_x)
        py = self.calc_position(node.y, self.min_y)

        # 检查当前位置是否在环境内
        if px < self.min_x or px > self.max_x:
            return False
        if py < self.min_x or py > self.max_y:
            return False

        return True

    @staticmethod
    def get_motion_model(x,y,cosfun):
        # dx, dy, cost
        data_x=x-1
        data_y=y-1
        model = [
            [0, 1, cosfun[data_x,data_y+1]],      # 上
            [0, -1,cosfun[data_x,data_y-1]],     # 下
            [-1, 0,cosfun[data_x-1,data_y]],     # 左
            [1, 0, cosfun[data_x+1,data_y]],      # 右
        ]
        return model

def ReadCsvData(filepath,xrow,ycol):
    file = open(filepath, "rb")
    filedata=np.loadtxt(file, delimiter=";")
    file.close()
    filearray=np.array(filedata)

    dataname = np.zeros(shape=(xrow,ycol))
    for i in range(xrow):
        for j in range(ycol):
            dataname[i,j]=filearray[i, j]
    return dataname
def main():
    datanamuber=50
    filepath = "E:/User/Desktop/datalab/u-wind1.csv"
    cosfun = ReadCsvData(filepath,datanamuber+2,datanamuber+2)
    # 设置起点，终点
    sx, sy = 3, 3
    gx, gy = 48, 48

    # 设置环境地图
    ox, oy = [], []

    # 设置四条边，
    for i in range(0, datanamuber):         # 下边
        ox.append(i)
        oy.append(0)
    for i in range(0, datanamuber):         # 左边
        ox.append(0)
        oy.append(i)
    for i in range(0, datanamuber):         # 上边
        ox.append(i)
        oy.append(datanamuber)
    for i in range(0, datanamuber):         # 右边
        ox.append(datanamuber)
        oy.append(i)

    if show:
        plt.plot(ox, oy, '.k')
        plt.plot(sx, sy, 'og')
        plt.plot(gx, gy, 'or')
        # plt.grid('True')
        plt.axis('equal')
        # plt.show()
    dijkstra = Dijkstra(ox, oy,cosfun,sx, sy)
    rx, ry = dijkstra.planning(sx, sy, gx, gy,cosfun)
    # for i in range(len(rx)):
    #     print(cosfun[rx[i]-1,ry[i]-1])
    # print(cosfun)
    if show:
        plt.plot(rx, ry, '-r')
        plt.pause(0.01)
        plt.show()

if __name__ == '__main__':
    show = True
    main()
