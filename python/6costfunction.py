import math
import numpy as np
import pandas as pd
## Function for importing data from Csv file
## numpy -> delimiter = ";"
def ReadCsvWind(filepath,xrow,ycol):
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
## pandas -> delimiter = "rows & columns"
def ReadCsvWave(filepath,xrow,ycol):
    file =pd.read_csv(filepath)

    num_row=xrow+2                  # +2: just in case for that: data is not enough for grid
    num_col=ycol+2
    dataname = file.values[0:num_row,0:num_col]
    return dataname
## Function for calculating wind resistance:
## R=0.5*rho_air*Avx*[Cair(relative wind angle)*V(relative wind velocity)^2-Cair(0)*V(ship velocity)^2]
def WindResCosFun(Vs,Axv,rho_air,wind_u,wind_v,Cair_extend):
    angle_ship=np.mat([0,90,180,270])  #Angle of ship: North, East, South and West
    for k in range(angle_ship.shape[1]):
        vel_wind_squre=np.mat(np.zeros(wind_u.shape))
        angle_wind=np.mat(np.zeros(wind_v.shape))
        angle_rel=np.mat(np.zeros(wind_u.shape))
        Cx=np.mat(np.zeros(wind_v.shape))
        if k == 0:
            wind_u_north=wind_u+Vs
            for i in range(wind_u.shape[0]):
                for j in range(wind_v.shape[1]):
                    vel_wind_squre[i,j]=wind_u_north[i, j]**2+wind_v[i, j]**2
                    if vel_wind_squre[i,j]==0:
                        angle_wind[i,j]=0
                    else:
                        angle_wind[i,j]=math.degrees(math.asin(wind_v[i,j]/math.sqrt(vel_wind_squre[i,j])))
                    
                    angle_rel[i,j]=angle_wind[i,j]-angle_ship[0,k]+180
                    while angle_rel[i,j]<0 or angle_rel[i,j]>=360:
                        if angle_rel[i,j]< 0:
                            angle_rel[i,j]=angle_rel[i,j]+360
                        if angle_rel[i,j]>=360:
                            angle_rel[i,j]=angle_rel[i,j]-360

                    weight=angle_rel[i,j]/10
                    index=int(angle_rel[i,j]//10)
                    Cx[i,j]=(weight-index)*Cair_extend[0,index+1]+(1-weight+index)*Cair_extend[0,index]
            Rwind_N=0.5*rho_air*Axv*(Cx*vel_wind_squre-Cair_extend[0,0]*Vs**2)
            CF_Wind_North=Rwind_N*Vs
        
        if k == 1:
            wind_v_east=wind_v+Vs
            for i in range(wind_u.shape[0]):
                for j in range(wind_v.shape[1]):
                    vel_wind_squre[i,j]=wind_u[i, j]**2+wind_v_east[i, j]**2
                    if vel_wind_squre[i,j]==0:
                        angle_wind[i,j]=0
                    else:
                        angle_wind[i,j]=math.degrees(math.asin(wind_v_east[i,j]/math.sqrt(vel_wind_squre[i,j])))
                    
                    angle_rel[i,j]=angle_wind[i,j]-angle_ship[0,k]+180
                    while angle_rel[i,j]<0 or angle_rel[i,j]>=360:
                        if angle_rel[i,j]< 0:
                            angle_rel[i,j]=angle_rel[i,j]+360
                        if angle_rel[i,j]>=360:
                            angle_rel[i,j]=angle_rel[i,j]-360

                    weight=angle_rel[i,j]/10
                    index=int(angle_rel[i,j]//10)
                    Cx[i,j]=(weight-index)*Cair_extend[0,index+1]+(1-weight+index)*Cair_extend[0,index]
            Rwind_E=0.5*rho_air*Axv*(Cx*vel_wind_squre-Cair_extend[0,0]*Vs**2)
            CF_Wind_East=Rwind_E*Vs

        if k == 2:
            wind_u_south=wind_u-Vs
            for i in range(wind_u.shape[0]):
                for j in range(wind_v.shape[1]):
                    vel_wind_squre[i,j]=wind_u_south[i,j]**2+wind_v[i,j]**2
                    if vel_wind_squre[i,j]==0:
                        angle_wind[i,j]=0
                    else:
                        angle_wind[i,j]=math.degrees(math.asin(wind_v[i,j]/math.sqrt(vel_wind_squre[i,j])))
                    
                    angle_rel[i,j]=angle_wind[i,j]-angle_ship[0,k]+180
                    while angle_rel[i,j]<0 or angle_rel[i,j]>=360:
                        if angle_rel[i,j]< 0:
                            angle_rel[i,j]=angle_rel[i,j]+360
                        if angle_rel[i,j]>=360:
                            angle_rel[i,j]=angle_rel[i,j]-360

                    weight=angle_rel[i,j]/10
                    index=int(angle_rel[i,j]//10)
                    Cx[i,j]=(weight-index)*Cair_extend[0,index+1]+(1-weight+index)*Cair_extend[0,index]
            Rwind_S=0.5*rho_air*Axv*(Cx*vel_wind_squre-Cair_extend[0,0]*Vs**2)
            CF_Wind_South=Rwind_S*Vs
        
        if k == 3:
            wind_v_west=wind_v-Vs
            for i in range(wind_u.shape[0]):
                for j in range(wind_v.shape[1]):
                    vel_wind_squre[i,j]=wind_u[i, j]**2+wind_v_west[i, j]**2
                    if vel_wind_squre[i,j]==0:
                        angle_wind[i,j]=0
                    else:
                        angle_wind[i,j]=math.degrees(math.asin(wind_v_west[i,j]/math.sqrt(vel_wind_squre[i,j])))
                    
                    angle_rel[i,j]=angle_wind[i,j]-angle_ship[0,k]+180
                    while angle_rel[i,j]<0 or angle_rel[i,j]>=360:
                        if angle_rel[i,j]< 0:
                            angle_rel[i,j]=angle_rel[i,j]+360
                        if angle_rel[i,j]>=360:
                            angle_rel[i,j]=angle_rel[i,j]-360

                    weight=angle_rel[i,j]/10
                    index=int(angle_rel[i,j]//10)
                    Cx[i,j]=(weight-index)*Cair_extend[0,index+1]+(1-weight+index)*Cair_extend[0,index]
            Rwind_W=0.5*rho_air*Axv*(Cx*vel_wind_squre-Cair_extend[0,0]*Vs**2)
            CF_Wind_West=Rwind_W*Vs

    return CF_Wind_North,CF_Wind_East,CF_Wind_South,CF_Wind_West
## Function for calculating total resistance (only constant velocity)
## Hollenbach Method (only return mean resistance)
def TotalResCosFun(Vs,L,Lwl,Los,B,TF,TA,CB,S,Dp,NRud,NBrac,NBoss,NThr,rho_sea,nu_sea,g):
   T = (TF+TA)/2

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

   Fn = Vs/((g*Lfn)**0.5) 	                    #Froude's number
   Fnkrit_help0=np.mat([1,CB,CB**2])          # Build Matrix for using transpose: Fnkrit_help0.T
   Fnkrit_help1 = d*Fnkrit_help0.T               # Fnkrit_help1=[[x]]   Matrix type
   Fnkrit=Fnkrit_help1[0,0]                           # Fnkrit=x                  Float type
   c1 = Fn/Fnkrit
   Rns = Vs*L/nu_sea						                   #Reynold's number for ship
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
   R_T_mean = C_Ts*rho_sea/2*Vs**2*S		   #Total resistance to the ship
   P_T_mean = R_T_mean*Vs	            	   #Propulsion power[W]

   return P_T_mean
## Function for calculating added resistance (-45 degrees ~ 45 degrees)
## R=1/16*rho_sea*g*H^2*B*sqrt(B/Lbwl)
def WaveResCosFun(Vs,Lbwl,B,rho_sea,wave_d,wave_h):
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
                P_wave=R_wave*Vs
        if k==0:
            CF_Wave_N=P_wave
        if k==1:
            CF_Wave_E=P_wave
        if k==2:
            CF_Wave_S=P_wave
        if k==0:
            CF_Wave_W=P_wave
    return CF_Wave_N,CF_Wave_E,CF_Wave_S,CF_Wave_W

## input some data
## for grids
## xrow,ycol
# xrow=int(input('number of rows:'))
# ycol=int(input('number of columns:'))
## for wind resistance
## Vs,rho_air,Cair,filepath1,filepath2
# Vs=float(input('Velocity of Ships(m/s):')) 
# rho_air=float(input('Density of Air (kg/m^3):')) 
# Cair=int(input('wind resistance coefficient(Be careful with Cair_extend):'))  
# filepath1=input('filepath of wind_u:')
# filepath2=input('filepath of wind_v:')
## for total resistance
## Vs,L,Lwl,Los,B,TF,TA,CB,S,Dp,NRud,NBrac,NBoss,NThr,rho_sea,nu_sea,g
# L=float(input('Length of Ship(m):'))
# Lwl=float(input('Length of Water Line(m):'))
# Los=float(input('Length over Surface(m):'))
# B=float(input('Beam(m):'))
# TF=float(input('Draft of Fore Propeller(m):'))
# TA=float(input('Draft of Aft Propeller(m):'))
# CB=float(input('Block coefficient:'))
# S=float(input('Wetted Surface(square m):'))
# Dp=float(input('Propeller diameter(m):'))
# NRud=float(input('Number of rudders:'))
# NBrac=float(input('Number of brackets:'))
# NBoss=float(input('Number of bossings:'))
# NThr=float(input('Number of side thrusters:'))
# rho_sea=float(input('Density of Sea Water (kg/m^3):'))
# nu_sea=float(input('Viscosity of Sea water (m/s^2):'))
# g=float(input('Gravitational Constant (m/s^2):'))
## for added resistance
## Vs,Lbwl,B,rho_sea,g,filepath3,filepath4
# Lbwl=float(input('Length of the Bow on the Water Line to 95% of maximum Beam (m):')
# filepath3=input('filepath of wave_d:')
# filepath4=input('filepath of wave_h:')

#data for test
xrow=1
ycol=1

## Wind Resistance
Vs=10                                       #Velocity of ship
Axv=1                                                #Area of maximum transverse section exposed to the wind

rho_air=1.293                                    #Density of sea water [kg/m^3]
filepath1 = "E:/User/Desktop/datalab/u-wind1.csv"
filepath2 = "E:/User/Desktop/datalab/v-wind1.csv"
wind_u = ReadCsvWind(filepath1,xrow,ycol)
wind_v = ReadCsvWind(filepath2,xrow,ycol)

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
    Cair_extend[0,i]=Cair[0,i]
    Cair_extend[0,2*Cair.size-i-2]=Cair[0,i]

CF_Wind_N,CF_Wind_E,CF_Wind_S,CF_Wind_W = WindResCosFun(Vs,Axv,rho_air,wind_u,wind_v,Cair_extend)

## Viscous/Friction+ Wave Resistance
L=6
Lwl=3
Los=5
B=2
TF=1.2
TA=1.25
CB=0.8
S=10
Dp=1
NRud=2
NBrac=1
NBoss=1
NThr=1

rho_sea = 1025
nu_sea = 1.1395E-6
g = 9.81

CF_m=TotalResCosFun(Vs,L,Lwl,Los,B,TF,TA,CB,S,Dp,NRud,NBrac,NBoss,NThr,rho_sea,nu_sea,g)

## Added Resistance
Lbwl=6      #length of the bow on the water line to 95% of maximum beam

filepath3 = "E:/User/Desktop/datalab/wave_mwd.csv"
filepath4 = "E:/User/Desktop/datalab/wave_swh.csv"
wave_d = ReadCsvWave(filepath3,xrow,ycol)
wave_h = ReadCsvWave(filepath4,xrow,ycol)

CF_Wave_N,CF_Wave_E,CF_Wave_S,CF_Wave_W=WaveResCosFun(Vs,Lbwl,B,rho_sea,wave_d,wave_h)

CF_N=CF_Wind_N+CF_Wave_N+CF_m
CF_E=CF_Wind_E+CF_Wave_E+CF_m
CF_S=CF_Wind_S+CF_Wave_S+CF_m
CF_W=CF_Wind_W+CF_Wave_W+CF_m