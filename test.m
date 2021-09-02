%% test
Vsvec=2.5;
L=6;
Lwl=3;
Los=5;
B=2;
TF=1.2;
TA=1.25;
CB=0.8;
S=10;
Dp=1;
NRud=2;
NBrac=1;
NBoss=1;
NThr=1;

T = (TF+TA)/2;
Fi = (CB/L) * (  (B/2) * (TF+TA)  )^.5;
k = .6 * Fi + 145 * Fi^3.5;

rho = 1025;
gravk = 9.81;   %Gravity
nu = 1.1395E-6; % viscosity


%Calculation of 'Froude length', Lfn:
if Los/L < 1
   Lfn = Los;
elseif (Los/L >= 1) & (Los/L < 1.1)
   Lfn = L+2/3*(Los-L);
elseif Los/L >= 1.1
   Lfn = 1.0667*L;   
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Constants from Hollenbachs paper: 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 'Mean' resistance coefficients
   a = [-0.3382 0.8086 -6.0258 -3.5632 9.4405 0.0146 0 0 0 0];  	%a1 means a(1) and so on
   b = [-0.57424 	 13.3893	90.5960; 	%b12 means b(1,2)
	   4.6614	-39.721	-351.483;
      -   1.14215	-12.3296	459.254];
   d = [0.854 -1.228 0.497];
   e = [2.1701 -0.1602];
   f = [0.17 0.20 0.60];
   g = [0.642 -0.635 0.150];

   % 'Minimum' resistance coefficients
   a_min = [-0.3382 0.8086 -6.0258 -3.5632 0 0 0 0 0 0];
   b_min = [-0.91424 13.3893 90.5960;...
         4.6614 -39.721 -351.483;...
         -1.14215 -12.3296 459.254];
   d_min = [0 0 0];
   e_min = [1 0];
   f_min = [0.17 0.2 0.6];
   g_min = [0.614 -0.717 0.261];


cc = 0;
% Loop over velocities
Vs = Vsvec;
	
	cc = cc + 1;
	
	% Froude's number
	Fn = Vs/sqrt(gravk*Lfn);			
	
	
	Fnkrit = d*[1 CB CB^2]';
	c1 = Fn/Fnkrit;
	c1_min = Fn/Fnkrit;
	
	Rns = Vs*L/nu;						% Reynold's number for ship
	CFs = 0.075/(log10(Rns)-2)^2;			% ITTC friction line for ship				

	% Calculation of C_R for given ship 
	% Mean value
	
	CRFnkrit = max(1.0,(Fn/Fnkrit)^c1);
	
	kL = e(1)*L^(e(2));
	
	% There is an error in the hollenbach paper and in Minsaas' 2003 textbook, which
	% is corrected in this formula by dividing by 10
	CRstandard = [1 CB CB^2]*(b*[1 Fn Fn^2]')/10;
	
	CR_hollenbach = CRstandard*CRFnkrit*kL*prod([T/B B/L Los/Lwl Lwl/L (1+(TA-TF)/L) ...
			Dp/TA (1+NRud) (1+NBrac) (1+NBoss) (1+NThr)].^a);
	
	CR = CR_hollenbach*B*T/S;   			% Resistance coefficient, scaled for wetted surface
	C_Ts = CFs + CR;				% Total resistance coeff. ship 
	R_T_mean = C_Ts*rho/2*Vs^2*S;			% Total resistance to the ship
	
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% Minimum values

	% There is an error in the hollenbach paper and in Minsaas' 2003 textbook, which
	% is corrected in this formula by dividing by 10
	CRstandard_min = [1 CB CB^2]*(b_min*[1 Fn Fn^2]')/10;
	test1=[T/B B/L Los/Lwl Lwl/L (1+(TA-TF)/L) Dp/TA (1+NRud) (1+NBrac) (1+NBoss) (1+NThr)];
    test2=test1.^a_min;
    test3=prod(test2);
	CR_hollenbach_min = CRstandard_min*test3;
	
	CR_min = CR_hollenbach_min*B*T/S;