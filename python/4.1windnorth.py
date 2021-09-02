import math
import numpy as np
#Function for importing data from Csv file
def ReadCsvData(filepath,xrow,ycol):
    file = open(filepath, "rb")
    filedata=np.loadtxt(file, delimiter=";")
    file.close()
    filearray=np.array(filedata)

    num_row=xrow+2                  # +2: just in case for that: data is not enough for grid
    num_col=ycol+2
    dataname = np.zeros(shape=(num_row,num_col))
    for i in range(num_row):
        for j in range(num_col):
            dataname[i,j]=filearray[i, j]
    return dataname
rho_air = 1.293             #Density of sea water [kg/m^3]
Axv=1                           #Area of maximum transverse section exposed to the wind
Cair=np.mat([-0.60, -0.87, -1.00, -1.00, -0.88, -0.85, -0.65, -0.42, -0.27, -0.09, 0.09, 0.49, 0.84, 1.39, 1.47, 1.34, 0.92, 0.82])  #Wind resistance coefficient: Data about General Cargo form ITTC
Cair_extend=np.mat(np.zeros((Cair.size*2-1)))        #Value Range of Inital Data:0-180 degrees 
Cair_extend[0,Cair.size-1]=Cair[0,Cair.size-1]          #Value Range of Extended Data:0-360 degrees 
for i in range(Cair.size-1):
    Cair_extend[0,i]=Cair[0,i]
    Cair_extend[0,2*Cair.size-i-2]=Cair[0,i]

xrow=1
ycol=1
filepath1 = "E:/User/Desktop/datalab/u-wind1.csv"
filepath2 = "E:/User/Desktop/datalab/v-wind1.csv"
wind_u = ReadCsvData(filepath1,xrow,ycol)
wind_v = ReadCsvData(filepath2,xrow,ycol)
vel_ship=10                                       #Velocity of ship
angle_ship=np.mat([0,90,180,270])  #Angle of ship: North, East, South and West
for k in range(angle_ship.shape[1]):
    vel_wind_squre=np.mat(np.zeros(wind_u.shape))
    angle_wind=np.mat(np.zeros(wind_v.shape))
    angle_rel=np.mat(np.zeros(wind_u.shape))
    Cx=np.mat(np.zeros(wind_v.shape))
    if k == 0:
        wind_u_north=wind_u+vel_ship
        for i in range(wind_u.shape[0]):
            for j in range(wind_v.shape[1]):
                vel_wind_squre[i,j]=wind_u_north[i, j]**2+wind_v[i, j]**2
                if vel_wind_squre[i,j]==0:
                    angle_wind[i,j]=0
                else:
                    angle_wind[i,j]=math.degrees(math.acos(wind_v[i,j]/math.sqrt(vel_wind_squre[i,j])))
                
                angle_rel[i,j]=angle_wind[i,j]-angle_ship[0,k]+180
                while angle_rel[i,j]<0 or angle_rel[i,j]>=360:
                    if angle_rel[i,j]< 0:
                        angle_rel[i,j]=angle_rel[i,j]+360
                    if angle_rel[i,j]>=360:
                        angle_rel[i,j]=angle_rel[i,j]-360
                print(math.sin(angle_rel[i,j]))
                weight=angle_rel[i,j]/10
                index=int(angle_rel[i,j]//10)
                Cx[i,j]=(weight-index)*Cair_extend[0,index+1]+(1-weight+index)*Cair_extend[0,index]
        Rwind_N=0.5*rho_air*Axv*(Cx*vel_wind_squre-Cair_extend[0,0]*vel_ship**2)
    # if k == 1:
    #     wind_v_east=wind_v+vel_ship
    #     for i in range(wind_u.shape[0]):
    #         for j in range(wind_v.shape[1]):
    #             vel_wind_squre[i,j]=wind_u[i, j]**2+wind_v_east[i, j]**2
    #             if vel_wind_squre[i,j]==0:
    #                 angle_wind[i,j]=0
    #             else:
    #                 angle_wind[i,j]=math.degrees(math.acos(wind_v_east[i,j]/vel_wind_squre[i,j]))
                
    #             angle_rel[i,j]=angle_wind[i,j]-angle_ship[0,k]+180
    #             while angle_rel[i,j]<0 or angle_rel[i,j]>=360:
    #                 if angle_rel[i,j]< 0:
    #                     angle_rel[i,j]=angle_rel[i,j]+360
    #                 if angle_rel[i,j]>=360:
    #                     angle_rel[i,j]=angle_rel[i,j]-360

    #             weight=angle_rel[i,j]/10
    #             index=int(angle_rel[i,j]//10)
    #             Cx[i,j]=(weight-index)*Cair_extend[0,index+1]+(1-weight+index)*Cair_extend[0,index]
    #     print(Cx)        
    # if k == 1:
    #     wind_v_east=wind_v+vel_ship
    #     for i in range(wind_u.shape[0]):
    #         for j in range(wind_v.shape[1]):
    #             vel_wind_squre[i,j]=wind_u[i, j]**2+wind_v_east[i, j]**2
    #             if vel_wind_squre[i,j]==0:
    #                 angle_wind[i,j]=0
    #             else:
    #                 angle_wind[i,j]=math.degrees(math.acos(wind_v_east[i,j]/vel_wind_squre[i,j]))
                
    #             angle_rel[i,j]=angle_wind[i,j]-angle_ship[0,k]+180
    #             while angle_rel[i,j]<0 or angle_rel[i,j]>=360:
    #                 if angle_rel[i,j]< 0:
    #                     angle_rel[i,j]=angle_rel[i,j]+360
    #                 if angle_rel[i,j]>=360:
    #                     angle_rel[i,j]=angle_rel[i,j]-360

    #             weight=angle_rel[i,j]/10
    #             index=int(angle_rel[i,j]//10)
    #             Cx[i,j]=(weight-index)*Cair_extend[0,index+1]+(1-weight+index)*Cair_extend[0,index]
print(Rwind_N)