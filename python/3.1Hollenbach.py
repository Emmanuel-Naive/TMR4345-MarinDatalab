import math
import numpy as np
#Multi-line comment shortcut: shirt+art+a
#test data for convenience
Vsvec=np.mat([1,0,2.5]) 
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

T = (TF+TA)/2
Fi = (CB/L) * (  (B/2) * (TF+TA)  )**0.5
k = 0.6 * Fi + 145 * Fi**3.5

rho = 1025         #Density of sea water [kg/m^3]
gravk = 9.81       #Gravitational constant [m/s^2]
nu = 1.1395E-6  #Viscosity of sea water [m/s^2]
""" 
Given in class
Rewritten by Weijian Yang, email: weijiany@stud.ntnu.no
Function: calculate cost caused by the reaction between hull and calm water
Method: Hollenbach
 """
#input some data
#Vsvec,L,Lwl,Los,B,TF,TA,CB,S,Dp,NRud,NBrac,NBoss,NThr
Vsvec=float(input('Velocity of Ships(m/s):'))               #The ship's velocity in m/s
L=float(input('Length of Ship(meters):'))                          #The ship's length between perpendiculars in meters
Lwl=float(input('Length of Water Line(meters):'))            #The ship's length in waterline in meters
Los=float(input('Length over Surface(meters):'))             #The ship's length over surface in meters
B=float(input('Beam(meters):'))                                        #The ship's breadth in meters
TF=float(input('Draft of Fore Propeller(meters):'))           #The ship's draught at front perpendiculars in meters
TA=float(input('Draft of Aft Propeller(meters):'))             #The ship's draught at after perpendicular in meters
CB=float(input('Block coefficient:'))                                  #The ship's block coefficient (dimensionless)
S=float(input('Wetted Surface(square meters):'))             #The ship's wetted surface area in square meters
Dp=float(input('Propeller diameter(meters):'))                 #Propeller diameter in meters
NRud=float(input('Number of rudders:'))                         #Number of rudders      
NBrac=float(input('Number of brackets:'))                       #Number of brackets
NBoss=float(input('Number of bossings:'))                      #Number of bossings
NThr=float(input('Number of side thrusters:'))                #Number of side thrusters


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
#'Minimum' resistance coefficients
a_min = np.mat([-0.3382,0.8086,-6.0258,-3.5632,0,0,0,0,0,0])
b_min = np.mat([[-0.91424,13.3893,90.5960],[4.6614,-39.721,-351.483],[-1.14215,-12.3296,459.254]])
d_min = np.mat([0,0,0])
e_min = np.mat([1,0])
f_min  = np.mat([0.17,0.2,0.6])
g_min = np.mat([0.614,-0.717,0.261])

cc = 0
# Loop over velocities
for i in range(Vsvec.size):
   Vs = Vsvec[0,i]       #Unit of Vs should be [m/s]

   cc = cc + 1

   Fn = Vs/((gravk*Lfn)**0.5) 	      #Froude's number
   Fnkrit_help0=np.mat([1,CB,CB**2])          # Build Matrix for using transpose: Fnkrit_help0.T
   Fnkrit_help1 = d*Fnkrit_help0.T               # Fnkrit_help1=[[x]]   Matrix type
   Fnkrit=Fnkrit_help1[0,0]                           # Fnkrit=x                  Float type
   c1 = Fn/Fnkrit
   c1_min = Fn/Fnkrit
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
   C_Ts = CFs + CR                                 #Total resistance coeff. ship 
   R_T_mean = C_Ts*rho/2*Vs**2*S		   #Total resistance to the ship

   #Minimum values

   #There is an error in the hollenbach paper and in Minsaas' 2003 textbook, which is corrected in this formula by dividing by 10
   CRstandard_min_help = Fnkrit_help0*(b_min*CRstandard_help0.T)/10
   CRstandard_min=CRstandard_min_help[0,0]
   #prod([T/B B/L Los/Lwl Lwl/L (1+(TA-TF)/L) Dp/TA (1+NRud) (1+NBrac) (1+NBoss) (1+NThr)].^a_min)
   prod_help_min1=np.mat(np.ones((1,10))) 
   for j in range(a_min.size):                                            #prod_help=[prod_help[0,i].^a_min[0,i]]
            prod_help_min1[0,j]=prod_help[0,j]**a_min[0,j]
   prod_help_min2=np.prod(prod_help_min1,axis = 1)                 #prod function
   prod_min=prod_help_min2[0,0]

   CR_hollenbach_min = CRstandard_min*prod_min
   CR_min = CR_hollenbach_min*B*T/S


   # Total resistance coefficient of the ship 
   C_Ts_min = CFs + CR_min
   # Total resistance	
   R_T_min = C_Ts_min*rho/2*Vs**2*S
   #Propulsion power
   P_E_mean = R_T_mean * Vs       #[W]
   P_E_min = R_T_min * Vs             #[W]

   #print sth.
   # print('Vs =', Vs/0.5144,'knots')
   # print('Mean values')
   # print('CRh:',CR)
   # print('CF:',CFs)
   # print('CT:',C_Ts)
   # print('RT:',R_T_mean,'N')
   # print('Minimum values')
   # print('CRh:',CR_min)
   # print('CF:',CFs)
   # print('CT:',C_Ts_min)
   # print('RT:',R_T_min)

   # % Store results for plotting
   CFsvec = np.mat(np.zeros((1,Vsvec.size)))
   CRvec = np.mat(np.zeros((1,Vsvec.size)))
   C_Tsvec = np.mat(np.zeros((1,Vsvec.size)))
   R_T_meanvec = np.mat(np.zeros((1,Vsvec.size)))
   CR_minvec = np.mat(np.zeros((1,Vsvec.size)))
   C_Ts_minvec = np.mat(np.zeros((1,Vsvec.size)))
   R_T_minvec = np.mat(np.zeros((1,Vsvec.size)))
   P_E_meanvec = np.mat(np.zeros((1,Vsvec.size)))
   P_E_minvec = np.mat(np.zeros((1,Vsvec.size)))
   CFsvec[0,i] = CFs
   CRvec[0,i] = CR
   C_Tsvec[0,i] = C_Ts
   R_T_meanvec[0,i] = R_T_mean
   CR_minvec[0,i] = CR_min
   C_Ts_minvec[0,i] = C_Ts_min
   R_T_minvec[0,i] = R_T_min
   P_E_meanvec[0,i] = P_E_mean
   P_E_minvec[0,i] = P_E_min

#This is the matrix the function returns
resistancedata =np.hstack((Vsvec.T,R_T_meanvec.T,R_T_minvec.T,P_E_meanvec.T,P_E_minvec.T))   
print(resistancedata)