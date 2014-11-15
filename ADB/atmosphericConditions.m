function [T, P, rhof] = atmosphericConditions(height)

%REFERENCE ATMOSPHERIC CONDITION

Temp = 320; % [K] Typical daily high in Green River in June
Pressure = 101325; %[Pa]
L = 0.0065; %Temperature lapse rate in [K/m]
R = 8.314; 
M = 0.0289644; %Molar mass of air in [kg/mol]
g = 9.81; 

%CALCULATE CONDITIONS

T = Temp - L * height; %Temperature at altitude [k]
P = Pressure * (1 - L*height/Temp)^(g*M/(R*L)); %Pressure at altitude [Pa]
rhof = P * M / (R * T); 