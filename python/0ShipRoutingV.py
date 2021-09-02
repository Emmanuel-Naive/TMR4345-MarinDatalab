import math
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
## Function for importing data from Csv file
## numpy -> delimiter = ";"
def ReadCsvWind(filepath,xcol,yrow):
    file = open(filepath, "rb")
    filedata=np.loadtxt(file, delimiter=";")
    file.close()
    filearray=np.array(filedata)
    # +2: just in case for that: data is not enough for grid
    num_row=xcol+2                  
    num_col=yrow+2
    dataname = np.zeros(shape=(num_row,num_col))
    for i in range(num_row):
        for j in range(num_col):
            dataname[i,j]=filearray[i, j]
    return dataname
## pandas -> delimiter = "rows & columns"
def ReadCsvWave(filepath,xcol,yrow):
    file =pd.read_csv(filepath)
    # +2: just in case for that: data is not enough for grid
    num_row=xcol+2 
    num_col=yrow+2
    dataname = file.values[0:num_row,0:num_col]
    return dataname
## Function for calculating wind resistance:
## R=0.5*rho_air*Avx*[Cair(relative wind angle)*V(relative wind velocity)^2-Cair(0)*V(ship velocity)^2]
def WindResCosFun(Vs,Axv,rho_air,wind_u,wind_v,Cair_extend,dt):
    angle_ship=np.mat([0,90,180,270])  #Angle of ship: North, East, South and West
    for k in range(angle_ship.shape[1]):
        vel_wind_squre=np.mat(np.zeros(wind_v.shape))
        vel_wt_squre=np.mat(np.zeros(wind_u.shape))
        angle_wt=np.mat(np.zeros(wind_u.shape))
        angle_wt_rel=np.mat(np.zeros(wind_v.shape))
        angle_rel=np.mat(np.zeros(wind_u.shape))
        Cx=np.mat(np.zeros(wind_v.shape))
        if k == 0:
            wind_u_north=wind_u+Vs
            for i in range(wind_u.shape[0]):
                for j in range(wind_v.shape[1]):
                    vel_wt_squre[i,j]=wind_u_north[i, j]**2+wind_v[i, j]**2     #Vwt:true wind velocity
                    if vel_wt_squre[i,j]==0:
                        angle_wt[i,j]=0
                    else:
                        angle_wt[i,j]=math.degrees(math.acos(wind_v[i,j]/math.sqrt(vel_wt_squre[i,j]))) #Awt:true wind direction
                    
                    angle_wt_rel[i,j]=angle_wt[i,j]-angle_ship[0,k]+180         #Awt-Aship
                    vel_wind_squre[i,j]=vel_wt_squre[i,j]+Vs**2+math.sqrt(vel_wt_squre[i,j])*Vs*math.cos(math.radians(angle_wt_rel[i,j])) #V_WRef
                    nume=math.sqrt(vel_wt_squre[i,j])*math.sin(math.radians(angle_wt_rel[i,j]))
                    deno=Vs+math.sqrt(vel_wt_squre[i,j])*math.cos(math.radians(angle_wt_rel[i,j]))
                    if deno==0:
                        if nume<0:
                            angle_rel[i,j]=-90
                        if nume>0:
                            angle_rel[i,j]=90
                        if nume==0:
                            angle_rel[i,j]=0
                    elif deno<0:
                        angle_rel[i,j]=math.degrees(math.atan(nume/deno))+180       #A_WRef
                    else:
                        angle_rel[i,j]=math.degrees(math.atan(nume/deno))
                    if angle_rel[i,j]< 0:
                        angle_rel[i,j]=angle_rel[i,j]+360
                    
                    weight=angle_rel[i,j]/10
                    index=int(angle_rel[i,j]//10)
                    Cx[i,j]=(weight-index)*Cair_extend[0,index+1]+(1-weight+index)*Cair_extend[0,index]
            Rwind_N=0.5*rho_air*Axv*(Cx*vel_wind_squre-Cair_extend[0,0]*Vs**2)
            CF_Wind_North=Rwind_N*Vs*dt

        if k == 1:
            wind_v_east=wind_v+Vs
            for i in range(wind_u.shape[0]):
                for j in range(wind_v.shape[1]):
                    vel_wt_squre[i,j]=wind_u[i, j]**2+wind_v_east[i, j]**2     #Vwt:true wind velocity
                    if vel_wt_squre[i,j]==0:
                        angle_wt[i,j]=0
                    else:
                        angle_wt[i,j]=math.degrees(math.asin(wind_v_east[i,j]/math.sqrt(vel_wt_squre[i,j]))) #Awt:true wind direction
                    
                    angle_wt_rel[i,j]=angle_wt[i,j]-angle_ship[0,k]+180         #Awt-Aship
                    vel_wind_squre[i,j]=vel_wt_squre[i,j]+Vs**2+math.sqrt(vel_wt_squre[i,j])*Vs*math.cos(math.radians(angle_wt_rel[i,j])) #V_WRef
                    nume=math.sqrt(vel_wt_squre[i,j])*math.sin(math.radians(angle_wt_rel[i,j]))
                    deno=Vs+math.sqrt(vel_wt_squre[i,j])*math.cos(math.radians(angle_wt_rel[i,j]))
                    if deno==0:
                        if nume<0:
                            angle_rel[i,j]=-90
                        if nume>0:
                            angle_rel[i,j]=90
                        if nume==0:
                            angle_rel[i,j]=0
                    elif deno<0:
                        angle_rel[i,j]=math.degrees(math.atan(nume/deno))+180       #A_WRef
                    else:
                        angle_rel[i,j]=math.degrees(math.atan(nume/deno))
                    if angle_rel[i,j]< 0:
                        angle_rel[i,j]=angle_rel[i,j]+360
                    
                    weight=angle_rel[i,j]/10
                    index=int(angle_rel[i,j]//10)
                    Cx[i,j]=(weight-index)*Cair_extend[0,index+1]+(1-weight+index)*Cair_extend[0,index]
            Rwind_E=0.5*rho_air*Axv*(Cx*vel_wind_squre-Cair_extend[0,0]*Vs**2)
            CF_Wind_East=Rwind_E*Vs*dt

        if k == 2:
            wind_u_south=wind_u-Vs
            for i in range(wind_u.shape[0]):
                for j in range(wind_v.shape[1]):
                    vel_wt_squre[i,j]=wind_u_south[i,j]**2+wind_v[i,j]**2     #Vwt:true wind velocity
                    if vel_wt_squre[i,j]==0:
                        angle_wt[i,j]=0
                    else:
                        angle_wt[i,j]=math.degrees(math.acos(wind_v[i,j]/math.sqrt(vel_wt_squre[i,j]))) #Awt:true wind direction
                    
                    angle_wt_rel[i,j]=angle_wt[i,j]-angle_ship[0,k]+180         #Awt-Aship
                    vel_wind_squre[i,j]=vel_wt_squre[i,j]+Vs**2+math.sqrt(vel_wt_squre[i,j])*Vs*math.cos(math.radians(angle_wt_rel[i,j])) #V_WRef
                    nume=math.sqrt(vel_wt_squre[i,j])*math.sin(math.radians(angle_wt_rel[i,j]))
                    deno=Vs+math.sqrt(vel_wt_squre[i,j])*math.cos(math.radians(angle_wt_rel[i,j]))
                    if deno==0:
                        if nume<0:
                            angle_rel[i,j]=-90
                        if nume>0:
                            angle_rel[i,j]=90
                        if nume==0:
                            angle_rel[i,j]=0
                    elif deno<0:
                        angle_rel[i,j]=math.degrees(math.atan(nume/deno))+180       #A_WRef
                    else:
                        angle_rel[i,j]=math.degrees(math.atan(nume/deno))
                    if angle_rel[i,j]< 0:
                        angle_rel[i,j]=angle_rel[i,j]+360
                    
                    weight=angle_rel[i,j]/10
                    index=int(angle_rel[i,j]//10)
                    Cx[i,j]=(weight-index)*Cair_extend[0,index+1]+(1-weight+index)*Cair_extend[0,index]
            Rwind_S=0.5*rho_air*Axv*(Cx*vel_wind_squre-Cair_extend[0,0]*Vs**2)
            CF_Wind_South=Rwind_S*Vs*dt
        
        if k == 3:
            wind_v_west=wind_v-Vs
            for i in range(wind_u.shape[0]):
                for j in range(wind_v.shape[1]):
                    vel_wt_squre[i,j]=wind_u[i, j]**2+wind_v_west[i, j]**2     #Vwt:true wind velocity
                    if vel_wt_squre[i,j]==0:
                        angle_wt[i,j]=0
                    else:
                        angle_wt[i,j]=math.degrees(math.asin(wind_v_west[i,j]/math.sqrt(vel_wt_squre[i,j]))) #Awt:true wind direction
                    
                    angle_wt_rel[i,j]=angle_wt[i,j]-angle_ship[0,k]+180         #Awt-Aship
                    vel_wind_squre[i,j]=vel_wt_squre[i,j]+Vs**2+math.sqrt(vel_wt_squre[i,j])*Vs*math.cos(math.radians(angle_wt_rel[i,j])) #V_WRef
                    nume=math.sqrt(vel_wt_squre[i,j])*math.sin(math.radians(angle_wt_rel[i,j]))
                    deno=Vs+math.sqrt(vel_wt_squre[i,j])*math.cos(math.radians(angle_wt_rel[i,j]))
                    if deno==0:
                        if nume<0:
                            angle_rel[i,j]=-90
                        if nume>0:
                            angle_rel[i,j]=90
                        if nume==0:
                            angle_rel[i,j]=0
                    elif deno<0:
                        angle_rel[i,j]=math.degrees(math.atan(nume/deno))+180       #A_WRef
                    else:
                        angle_rel[i,j]=math.degrees(math.atan(nume/deno))
                    if angle_rel[i,j]< 0:
                        angle_rel[i,j]=angle_rel[i,j]+360
                    
                    weight=angle_rel[i,j]/10
                    index=int(angle_rel[i,j]//10)
                    Cx[i,j]=(weight-index)*Cair_extend[0,index+1]+(1-weight+index)*Cair_extend[0,index]
            Rwind_W=0.5*rho_air*Axv*(Cx*vel_wind_squre-Cair_extend[0,0]*Vs**2)
            CF_Wind_West=Rwind_W*Vs*dt

    return CF_Wind_North,CF_Wind_East,CF_Wind_South,CF_Wind_West
## Function for calculating total resistance (only constant velocity)
## Hollenbach Method (only return mean resistance)
def TotalResCosFun(Vs,L,Lwl,Los,B,TF,TA,CB,S,Dp,NRud,NBrac,NBoss,NThr,rho_sea,nu_sea,g,dt):
   T = (TF+TA)/2

   rho = 1025         #Density of sea water [kg/m^3]
   gravk = 9.81       #Gravitational constant [m/s^2]
   nu = 1.1395E-6  #Viscosity of sea water [m/s^2]

   #Calculation of 'Froude length', Lfn:
   if Los/L < 1:
      Lfn = Los
   elif (Los/L >= 1) & (Los/L < 1.1):
      Lfn = L+2/3*(Los-L)
   elif Los/L >= 1.1:
      Lfn = 1.0667*L

   # 'Mean' resistance coefficients
   a = np.mat([-0.3382, 0.8086, -6.0258, -3.5632, 9.4405, 0.0146, 0, 0, 0, 0])                                                #a1 means a[0,0]
   b = np.mat([[-0.57424, 13.3893, 90.5960],[4.6614, -39.721, -351.483],[-1.14215, -12.3296, 459.254]]) 	#b12 means b[0,1]
   d = np.mat([0.854, -1.228, 0.497])
   e = np.mat([2.1701, -0.1602])
   f  = np.mat([0.17, 0.20, 0.60])
   g = np.mat([0.642, -0.635, 0.150])

   Fn = Vs/((gravk*Lfn)**0.5) 	                    #Froude's number
   Fnkrit_help0=np.mat([1,CB,CB**2])          # Build Matrix for using transpose: Fnkrit_help0.T
   Fnkrit_help1 = d*Fnkrit_help0.T               # Fnkrit_help1=[[x]]   Matrix type
   Fnkrit=Fnkrit_help1[0,0]                           # Fnkrit=x                  Float type
   c1 = Fn/Fnkrit
   Rns = Vs*L/nu						                          #Reynold's number for ship
   if Rns == 0:                                                      #Rns=0,log would get stuck
      CFs =0
   else :
      CFs = 0.075/(math.log10(Rns)-2)**2			#ITTC friction line for ship

   
   # Calculation of C_R for given ship 
   # Mean value
   CRFnkrit = max(1.0,(Fn/Fnkrit)**c1)
   kL = e[0,0]*L**(e[0,1])

   # There is an error in the hollenbach paper and in Minsaas' 2003 textbook, which is corrected in this formula by dividing by 10
   CRstandard_help0=np.mat([1,Fn,Fn**2])
   CRstandard_help1=Fnkrit_help0*(b*CRstandard_help0.T)/10
   CRstandard=CRstandard_help1[0,0]

   #prod([T/B B/L Los/Lwl Lwl/L (1+(TA-TF)/L) Dp/TA (1+NRud) (1+NBrac) (1+NBoss) (1+NThr)].^a)
   prod_help=np.mat([T/B,B/L,Los/Lwl,Lwl/L,1+(TA-TF)/L,Dp/TA,1+NRud,1+NBrac,1+NBoss,1+NThr])
   prod_help1=np.mat(np.ones((1,10)))                          #build a Matrix[1,10]
   for j in range(a.size):                                            #prod_help=[prod_help[0,i].^a_min[0,i]]
            prod_help1[0,j]=prod_help[0,j]**a[0,j]
   prod_help2=np.prod(prod_help1,axis = 1)                 #prod function
   prod=prod_help2[0,0]

   CR_hollenbach = CRstandard*CRFnkrit*kL*prod
   CR = CR_hollenbach*B*T/S  			       #Resistance coefficient, scaled for wetted surface
   C_Ts = CFs + CR                                     #Total resistance coeff. ship 
   R_T_mean = C_Ts*rho/2*Vs**2*S		   #Total resistance to the ship
   
   Fn_min=min(f[0,0],f[0,0]+f[0,1]*(f[0,2]-CB))
   Fn_max=g[0,0]+g[0,1]*CB+g[0,2]*CB**3

   if Fn>Fn_max:
      #R_T=h1*R_T_mean 
      R_T=1.204*R_T_mean  #h1=1.204
   elif Fn<Fn_min:
      #CRFnkrit=kL=1.0
      CR_min=CRstandard*prod*B*T/S
      R_T= (CFs + CR_min)*rho/2*Vs**2*S
   else:
      R_T=R_T_mean
   
   C_T = R_T*Vs*dt
   return C_T
## Function for calculating added resistance (-45 degrees ~ 45 degrees)
## R=1/16*rho_sea*g*H^2*B*sqrt(B/Lbwl)
def WaveResCosFun(Vs,Lbwl,B,rho_sea,wave_d,wave_h,dt):
    angle_ship=np.mat([0,90,180,270])  #Angle of ship: North, East, South and West
    for k in range(angle_ship.shape[1]):
        angle_rel=np.mat(np.zeros(wave_d.shape))
        R_wave=np.mat(np.zeros(wave_h.shape))
        for i in range(wave_d.shape[0]):
            for j in range(wave_h.shape[1]):
                angle_rel[i,j]=wave_d[i,j]-angle_ship[0,k]
                while angle_rel[i,j]<-180 or angle_rel[i,j]>=180:
                    if angle_rel[i,j]< -180:
                        angle_rel[i,j]=angle_rel[i,j]+360
                    if angle_rel[i,j]>=180:
                        angle_rel[i,j]=angle_rel[i,j]-360
                
                if angle_rel[i,j]<-45 or angle_rel[i,j]>45:
                    R_wave[i,j]=0
                else:
                    R_wave[i,j]=rho_sea*g*B*math.sqrt(B/Lbwl)/16*wave_h[i,j]**2
                C_wave=R_wave*Vs*dt
        if k==0:
            CF_Wave_N=C_wave
        if k==1:
            CF_Wave_E=C_wave
        if k==2:
            CF_Wave_S=C_wave
        if k==0:
            CF_Wave_W=C_wave
    return CF_Wave_N,CF_Wave_E,CF_Wave_S,CF_Wave_W
## Dijkstra's alogrithm (find path)

class Dijkstra:

    def __init__(self, ox, oy,sx,sy,CF_N,CF_E,CF_S,CF_W):
        # 初始化地图的情况
        self.min_x = None
        self.max_x = None
        self.min_y = None
        self.max_y = None
        self.x_grid_num = None
        self.y_grid_num = None
        self.obstacle_map = None

        self.calc_obstacle_grid_map(ox, oy)         # 构建环境栅格地图
        self.robot_motion = self.get_motion_model(sx,sy,CF_N,CF_E,CF_S,CF_W)

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

    def planning(self, sx, sy, gx, gy,CF_N,CF_E,CF_S,CF_W):
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
                goal_node.parent_index = current.parent_index
                goal_node.cost = current.cost
                break

            # (3). 将该节点从 open_set 中取出，并加入到 close_set 中
            del open_set[c_id]
            close_set[c_id] = current

            # (4). 根据机器人的运动模式，在栅格地图中探索当前位置出发到达的下一可能位置
            self.robot_motion = self.get_motion_model(current.x,current.y,CF_N,CF_E,CF_S,CF_W)
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

        return rx, ry, goal_node.cost

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
    def get_motion_model(x,y,CF_N,CF_E,CF_S,CF_W):
        # dx, dy, cost
        data_x=x
        data_y=y
        model = [
            [0, 1, CF_N[data_x,data_y+1]],         # North
            [0, -1,CF_S[data_x,data_y-1]],          # South
            [-1, 0,CF_E[data_x-1,data_y]],          # East
            [1, 0, CF_W[data_x+1,data_y]],        # West
        ]
        return model

# ## for path planning
# ## t,sx,sy,gx,gy,xcol,yrow,num_iter
# t=float(input('given time(s):'))
# sx=int(input('start point(x):'))
# sy=int(input('start point(y):'))
# gx=int(input('goal point(x):'))
# gy=int(input('goal point(y):'))
# # it is a better choice to choose xcol=yrow
# xcol=int(input('number of columns(>abs(sy-gy)):'))
# yrow=int(input('number of rows(>abs(sx-gx):'))
# num_iter=int(input('number of iterations:'))
# ## for wind resistance
# ## Vs,rho_air,Cair,filepath1,filepath2
# Axv=float(input('Area of maximum transverse section exposed to the wind(m^2):'))
# rho_air=float(input('Density of Air (kg/m^3):')) 
# Cair=np.mat(float(input('wind resistance coefficient(Be careful with Cair_extend):')))
# filepath1=input('filepath of wind_u:')
# filepath2=input('filepath of wind_v:')
# ## for total resistance
# ## Vs,L,Lwl,Los,B,TF,TA,CB,S,Dp,NRud,NBrac,NBoss,NThr,rho_sea,nu_sea,g
# L=float(input('Length of Ship(m):'))
# Lwl=float(input('Length of Water Line(m):'))
# Los=float(input('Length over Surface(m):'))
# B=float(input('Beam(m):'))
# TF=float(input('Draft of Fore Propeller(m):'))
# TA=float(input('Draft of Aft Propeller(m):'))
# CB=float(input('Block coefficient:'))
# S=float(input('Wetted Surface(m^2):'))
# Dp=float(input('Propeller diameter(m):'))
# NRud=float(input('Number of rudders:'))
# NBrac=float(input('Number of brackets:'))
# NBoss=float(input('Number of bossings:'))
# NThr=float(input('Number of side thrusters:'))
# rho_sea=float(input('Density of Sea Water (kg/m^3):'))
# nu_sea=float(input('Viscosity of Sea water (m/s^2):'))
# g=float(input('Gravitational Constant (m/s^2):'))
# ## for added resistance
# ## Vs,Lbwl,B,rho_sea,g,filepath3,filepath4
# Lbwl=float(input('Length of the Bow on the Water Line to 95% of maximum Beam (m):')"
# filepath3=input('filepath of wave_d:')
# filepath4=input('filepath of wave_h:')

#data for test
t=100                                    #Given Time [s]
sx, sy = 5, 3                        #x,y of Start Point
gx, gy = 13, 17                   #x,y of Goal Point
xcol=20                             #Number of Columns
yrow=20                            #Number of Rows

## Wind Resistance

Axv=1                                                #Area of Maximum Yransverse Section Exposed to the Wind [m^2]

rho_air=1.293                                    #Density of Sea Water [kg/m^3]
filepath1 = "E:/User/Desktop/datalab/u-wind1.csv"
filepath2 = "E:/User/Desktop/datalab/v-wind1.csv"
wind_u = ReadCsvWind(filepath1,xcol,yrow)
wind_v = ReadCsvWind(filepath2,xcol,yrow)

""" 
If a new Cair is used, be careful with below functions for Cair_extend
Inital data about General Cargo form ITTC (Value range: 0-180 degrees): 
Cair=[-0.60, -0.87, -1.00, -1.00, -0.88, -0.85, -0.65, -0.42, -0.27, -0.09, 0.09, 0.49, 0.84, 1.39, 1.47, 1.34, 0.92, 0.82]
Step 1: Wind resistance coefficient: Data about General Cargo form ITTC, but no data for 130 degree angle
Cair=np.mat([-0.60, -0.87, -1.00, -1.00, -0.88, -0.85, -0.65, -0.42, -0.27, -0.09, 0.09, 0.49, 0.84, Cair_130,1.39, 1.47, 1.34, 0.92, 0.82])
Step 2: Extend this matrix and each element would be changed as its own opposite value (Value range: 0-360 degrees):
"""
# Step 1
Cair_130=0.5*0.84+0.5*1.39                                     #calculate the coefficient of 130 degree angle with that of 120 and 140 degree angle
Cair=np.mat([-0.60, -0.87, -1.00, -1.00, -0.88, -0.85, -0.65, -0.42, -0.27, -0.09, 0.09, 0.49, 0.84, Cair_130,1.39, 1.47, 1.34, 0.92, 0.82])
# Step 2
Cair_extend=np.mat(np.zeros((Cair.size*2-1)))        #Value Range of Inital Data:0-180 degrees 
Cair_extend[0,Cair.size-1]=Cair[0,Cair.size-1]          #Value Range of Extended Data:0-360 degrees 
for i in range(Cair.size-1):
    Cair_extend[0,i]=-Cair[0,i]
    Cair_extend[0,2*Cair.size-i-2]=-Cair[0,i]

## Viscous/Friction+ Wave Resistance
L=6                              #Length of Ship(m)
Lwl=3                          #Length of Water Line(m)
Los=5                          #Length over Surface(m)
B=2                              #Beam(m)
TF=1.2                         #Draft of Fore Propeller(m)
TA=1.25                       #Draft of Aft Propeller(m)
CB=0.8                         #Block coefficient
S=10                            #Wetted Surface(m^2)
Dp=1                           #Propeller diameter(m)
NRud=2                       #Number of rudders
NBrac=1                       #Number of brackets
NBoss=1                       #Number of bossings
NThr=1                         #Number of side thrusters

rho_sea = 1025              #Density of Sea Water (kg/m^3)
nu_sea = 1.1395E-6       #Viscosity of Sea water (m/s^2)
g = 9.81                         #Gravitational Constant (m/s^2)

## Added Resistance
Lbwl=6                          #Length of the Bow on the Water Line to 95% of Maximum Beam

filepath3 = "E:/User/Desktop/datalab/wave_mwd.csv"
filepath4 = "E:/User/Desktop/datalab/wave_swh.csv"
wave_d = ReadCsvWave(filepath3,xcol,yrow)
wave_h = ReadCsvWave(filepath4,xcol,yrow)

# 设置环境地图
ox, oy = [], []

# 设置四条边，
for i in range(0, xcol):        #north boundary
    ox.append(i)
    oy.append(xcol)
for i in range(0, yrow):         #east boundary
    ox.append(0)
    oy.append(i)
for i in range(0, xcol):        #south boundary
    ox.append(i)
    oy.append(0)
for i in range(0, yrow):         #west boundary
    ox.append(yrow)
    oy.append(i)

# num_iter(max)=xcol*yrow
num_iter=30
T_cost=np.mat(np.zeros((1,num_iter)))
TotalCost=T_cost
""" 
find minimum cost in each velocity
Vs_min=(abs(sx-gx)+abs(sy-gy)/t
Vs_max=2*Vs_min
 """
show = False
for i in range(num_iter):
    length=(abs(sx-gx)+abs(sy-gy))+i
    Vs=float(length/t)
    delta_T=1/Vs
    CF_Wind_N,CF_Wind_E,CF_Wind_S,CF_Wind_W = WindResCosFun(Vs,Axv,rho_air,wind_u,wind_v,Cair_extend,delta_T)
    CF_m=TotalResCosFun(Vs,L,Lwl,Los,B,TF,TA,CB,S,Dp,NRud,NBrac,NBoss,NThr,rho_sea,nu_sea,g,delta_T)
    CF_Wave_N,CF_Wave_E,CF_Wave_S,CF_Wave_W=WaveResCosFun(Vs,Lbwl,B,rho_sea,wave_d,wave_h,delta_T)
    CF_N=CF_Wind_N+CF_Wave_N+CF_m
    CF_E=CF_Wind_E+CF_Wave_E+CF_m
    CF_S=CF_Wind_S+CF_Wave_S+CF_m
    CF_W=CF_Wind_W+CF_Wave_W+CF_m
    dijkstra = Dijkstra(ox, oy,sx, sy,CF_N,CF_E,CF_S,CF_W)
    rx, ry, T_cost[0,i] = dijkstra.planning(sx, sy, gx, gy,CF_N,CF_E,CF_S,CF_W)
    if len(rx)*delta_T>t:
        print('The ship cannot reach the goal point in given time with this speed',Vs ,' [m/s]')
        TotalCost[0,i]=0
    else:
        print('The ship could reach the goal point in given time with this speed',Vs ,' [m/s]')
        TotalCost[0,i]=T_cost[0,i]
if np.max(TotalCost)==0:
    print('The effective path could not be found after ',num_iter,' iterations' )
else:
    minx,miny= np.where(TotalCost == np.min(TotalCost[np.nonzero(TotalCost)]))
    length_E=(abs(sx-gx)+abs(sy-gy))+miny
    Vs_E=float(length_E/t)
    delta_TE=1/Vs_E
    CFE_Wind_N,CFE_Wind_E,CFE_Wind_S,CFE_Wind_W = WindResCosFun(Vs_E,Axv,rho_air,wind_u,wind_v,Cair_extend,delta_T)
    CFE_m=TotalResCosFun(Vs_E,L,Lwl,Los,B,TF,TA,CB,S,Dp,NRud,NBrac,NBoss,NThr,rho_sea,nu_sea,g,delta_T)
    CFE_Wave_N,CFE_Wave_E,CFE_Wave_S,CFE_Wave_W=WaveResCosFun(Vs_E,Lbwl,B,rho_sea,wave_d,wave_h,delta_T)
    CFE_N=CFE_Wind_N+CFE_Wave_N+CFE_m
    CFE_E=CFE_Wind_E+CFE_Wave_E+CFE_m
    CFE_S=CFE_Wind_S+CFE_Wave_S+CFE_m
    CFE_W=CFE_Wind_W+CFE_Wave_W+CFE_m
    show = True
    if show:
        plt.plot(ox, oy, '.k')
        plt.plot(sx, sy, 'og')
        plt.plot(gx, gy, 'or')
        # plt.grid('True')
        plt.axis('equal')
        # plt.show()
    dijkstra = Dijkstra(ox, oy,sx, sy,CFE_N,CFE_E,CFE_S,CFE_W)
    rx, ry, TotalCost_min = dijkstra.planning(sx, sy, gx, gy,CFE_N,CFE_E,CFE_S,CFE_W)
    if show:
        print('Find the effective path' )
        print('The efficient velocity is [m/s]:',Vs_E)
        print('The minimum cost is [J]:',TotalCost_min)
        plt.plot(rx, ry, '-r')
        plt.pause(0.01)
        plt.show()