function [T, P, rho, mu] = atmosphericProperties(height, T0, h0)

%%%%%%%%%%%%%%%%%%%%%% REFERENCE ATMOSPHERIC CONDITION %%%%%%%%%%%%%%%%%%%%

P_ref = 101325; 
L = 0.0065; %Temperature lapse rate in K/m
T_ref = T0 + L*h0; 
R = 8.314; 
M = 0.0289644; %Molar mass of air in kg/mol
g = 9.81; 

%% Sutherland's Constants %%

mu0_S = 18.27e-6; %Reference viscosity (Pa s)
C_S = 120; %Constant (K)
T0_S = 291.15; %Temperature (K)

%%%%%%%%%%%%%%%%%%%%%%%%% CALCULATE NEW CONDITIONS %%%%%%%%%%%%%%%%%%%%%%%%

T = T_ref - L * height; %Temperature at altitude
P = P_ref * (1 - L*height/T_ref)^(g*M/(R*L)); %Pressure at altitude
rho = P * M / (R * T); 
mu = mu0_S * (T0_S + C_S)/(T+C_S) * (T/T0_S)^(3/2); 

