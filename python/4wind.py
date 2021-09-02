import math
import numpy as np
## Function for importing data from Csv file
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
## Function for calculating wind resistance:
## R=0.5*rho_air*Avx*[Cair(relative wind angle)*V(relative wind velocity)^2-Cair(0)*V(ship velocity)^2]
def WindResCosFun(vel_ship,Axv,wind_u,wind_v,Cair_extend):
    rho_air=1.293                                    #Density of sea water [kg/m^3]
    angle_ship=np.mat([0,90,180,270])  #Angle of ship: North, East, South and West
    for k in range(angle_ship.shape[1]):
        vel_wind_squre=np.mat(np.zeros(wind_v.shape))
        vel_wt_squre=np.mat(np.zeros(wind_u.shape))
        angle_wt=np.mat(np.zeros(wind_u.shape))
        angle_wt_rel=np.mat(np.zeros(wind_v.shape))
        angle_rel=np.mat(np.zeros(wind_u.shape))
        Cx=np.mat(np.zeros(wind_v.shape))
        if k == 0:
            wind_u_north=wind_u+vel_ship
            for i in range(wind_u.shape[0]):
                for j in range(wind_v.shape[1]):
                    vel_wt_squre[i,j]=wind_u_north[i, j]**2+wind_v[i, j]**2     #Vwt:true wind velocity
                    if vel_wt_squre[i,j]==0:
                        angle_wt[i,j]=0
                    else:
                        angle_wt[i,j]=math.degrees(math.acos(wind_v[i,j]/math.sqrt(vel_wt_squre[i,j]))) #Awt:true wind direction
                    
                    angle_wt_rel[i,j]=angle_wt[i,j]-angle_ship[0,k]+180         #Awt-Aship
                    vel_wind_squre[i,j]=vel_wt_squre[i,j]+vel_ship**2+math.sqrt(vel_wt_squre[i,j])*vel_ship*math.cos(math.radians(angle_wt_rel[i,j])) #V_WRef
                    nume=math.sqrt(vel_wt_squre[i,j])*math.sin(math.radians(angle_wt_rel[i,j]))
                    deno=vel_ship+math.sqrt(vel_wt_squre[i,j])*math.cos(math.radians(angle_wt_rel[i,j]))
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
                    print(angle_rel[i,j])
                    weight=angle_rel[i,j]/10
                    index=int(angle_rel[i,j]//10)
                    Cx[i,j]=(weight-index)*Cair_extend[0,index+1]+(1-weight+index)*Cair_extend[0,index]
            Rwind_N=0.5*rho_air*Axv*(Cx*vel_wind_squre-Cair_extend[0,0]*vel_ship**2)
            CF_Wind_North=Rwind_N*vel_ship
        if k == 1:
            wind_v_east=wind_v+vel_ship
            for i in range(wind_u.shape[0]):
                for j in range(wind_v.shape[1]):
                    vel_wt_squre[i,j]=wind_u[i, j]**2+wind_v_east[i, j]**2     #Vwt:true wind velocity
                    if vel_wt_squre[i,j]==0:
                        angle_wt[i,j]=0
                    else:
                        angle_wt[i,j]=math.degrees(math.asin(wind_v_east[i,j]/math.sqrt(vel_wt_squre[i,j]))) #Awt:true wind direction
                    
                    angle_wt_rel[i,j]=angle_wt[i,j]-angle_ship[0,k]+180         #Awt-Aship
                    vel_wind_squre[i,j]=vel_wt_squre[i,j]+vel_ship**2+math.sqrt(vel_wt_squre[i,j])*vel_ship*math.cos(math.radians(angle_wt_rel[i,j])) #V_WRef
                    nume=math.sqrt(vel_wt_squre[i,j])*math.sin(math.radians(angle_wt_rel[i,j]))
                    deno=vel_ship+math.sqrt(vel_wt_squre[i,j])*math.cos(math.radians(angle_wt_rel[i,j]))
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
            Rwind_E=0.5*rho_air*Axv*(Cx*vel_wind_squre-Cair_extend[0,0]*vel_ship**2)
            CF_Wind_East=Rwind_E*vel_ship

        if k == 2:
            wind_u_south=wind_u-vel_ship
            for i in range(wind_u.shape[0]):
                for j in range(wind_v.shape[1]):
                    vel_wt_squre[i,j]=wind_u_south[i,j]**2+wind_v[i,j]**2     #Vwt:true wind velocity
                    if vel_wt_squre[i,j]==0:
                        angle_wt[i,j]=0
                    else:
                        angle_wt[i,j]=math.degrees(math.acos(wind_v[i,j]/math.sqrt(vel_wt_squre[i,j]))) #Awt:true wind direction
                    
                    angle_wt_rel[i,j]=angle_wt[i,j]-angle_ship[0,k]+180         #Awt-Aship
                    vel_wind_squre[i,j]=vel_wt_squre[i,j]+vel_ship**2+math.sqrt(vel_wt_squre[i,j])*vel_ship*math.cos(math.radians(angle_wt_rel[i,j])) #V_WRef
                    nume=math.sqrt(vel_wt_squre[i,j])*math.sin(math.radians(angle_wt_rel[i,j]))
                    deno=vel_ship+math.sqrt(vel_wt_squre[i,j])*math.cos(math.radians(angle_wt_rel[i,j]))
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
            Rwind_S=0.5*rho_air*Axv*(Cx*vel_wind_squre-Cair_extend[0,0]*vel_ship**2)
            CF_Wind_South=Rwind_S*vel_ship
        
        if k == 3:
            wind_v_west=wind_v-vel_ship
            for i in range(wind_u.shape[0]):
                for j in range(wind_v.shape[1]):
                    vel_wt_squre[i,j]=wind_u[i, j]**2+wind_v_west[i, j]**2     #Vwt:true wind velocity
                    if vel_wt_squre[i,j]==0:
                        angle_wt[i,j]=0
                    else:
                        angle_wt[i,j]=math.degrees(math.asin(wind_v_west[i,j]/math.sqrt(vel_wt_squre[i,j]))) #Awt:true wind direction
                    
                    angle_wt_rel[i,j]=angle_wt[i,j]-angle_ship[0,k]+180         #Awt-Aship
                    vel_wind_squre[i,j]=vel_wt_squre[i,j]+vel_ship**2+math.sqrt(vel_wt_squre[i,j])*vel_ship*math.cos(math.radians(angle_wt_rel[i,j])) #V_WRef
                    nume=math.sqrt(vel_wt_squre[i,j])*math.sin(math.radians(angle_wt_rel[i,j]))
                    deno=vel_ship+math.sqrt(vel_wt_squre[i,j])*math.cos(math.radians(angle_wt_rel[i,j]))
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
            Rwind_W=0.5*rho_air*Axv*(Cx*vel_wind_squre-Cair_extend[0,0]*vel_ship**2)
            CF_Wind_West=Rwind_W*vel_ship

    return CF_Wind_North,CF_Wind_East,CF_Wind_South,CF_Wind_West

## for grids
# xrow,ycol,filepath1,filepath2
# xrow=int(input('number of rows:'))
# ycol=int(input('number of columns:'))
# filepath1=input('filepath of wind_u:')
# filepath2=input('filepath of wind_v:')
## for wind resistance
# Cair,vel_ship
# Cair=int(input('wind resistance coefficient(Be careful with Cair_extend):'))                  
# vel_ship=float(input('Velocity of Ships(m/s):')) 

#data for test
xrow=1
ycol=1
filepath1 = "E:/User/Desktop/datalab/u-wind1.csv"
filepath2 = "E:/User/Desktop/datalab/v-wind1.csv"

vel_ship=15                                    #Velocity of ship
Axv=1                                                #Area of maximum transverse section exposed to the wind

wind_u = ReadCsvData(filepath1,xrow,ycol)
wind_v = ReadCsvData(filepath2,xrow,ycol)

""" 
If a new Cair is used, be careful with below functions for Cair_extend
Inital data about General Cargo form ITTC (Value range: 0-180 degrees): 
Cair=[-0.60, -0.87, -1.00, -1.00, -0.88, -0.85, -0.65, -0.42, -0.27, -0.09, 0.09, 0.49, 0.84, 1.39, 1.47, 1.34, 0.92, 0.82]
Step 1: Wind resistance coefficient: Data about General Cargo form ITTC, but no data for 130 degree angle
Cair=np.mat([-0.60, -0.87, -1.00, -1.00, -0.88, -0.85, -0.65, -0.42, -0.27, -0.09, 0.09, 0.49, 0.84, Cair_130,1.39, 1.47, 1.34, 0.92, 0.82])
Step 2: Extend this matrix (Value range: 0-360 degrees):
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

CF_W_N,CF_W_E,CF_W_S,CF_W_W = WindResCosFun(vel_ship,Axv,wind_u,wind_v,Cair_extend)
print(wind_u)
print(wind_v)

print(CF_W_N)
print(CF_W_S)