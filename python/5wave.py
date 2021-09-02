import math
import numpy as np
import pandas as pd
def ReadCsvWave(filepath,xrow,ycol):
    file =pd.read_csv(filepath)

    num_row=xrow+2                  # +2: just in case for that: data is not enough for grid
    num_col=ycol+2
    dataname = file.values[0:num_row,0:num_col]
    return dataname

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

xrow=1
ycol=1

Vs=10
Lbwl=6      #length of the bow on the water line to 95% of maximum beam
B=2

rho_sea = 1025
g = 9.81
filepath3 = "E:/User/Desktop/datalab/wave_mwd.csv"
filepath4 = "E:/User/Desktop/datalab/wave_swh.csv"
wave_d = ReadCsvWave(filepath3,xrow,ycol)
wave_h = ReadCsvWave(filepath4,xrow,ycol)

#R=1/16*rho*g*H^2*B*sqrt(B/L)

CF_Wave_N,CF_Wave_E,CF_Wave_S,CF_Wave_W=WaveResCosFun(Vs,Lbwl,B,rho_sea,wave_d,wave_h)