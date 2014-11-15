function [Cd D M Re] = Aerodynamics(v, T, mu, rho, D_o)

c = importdata('DragCoefficients.txt'); 
M_ref = c.data(:,1); Cd_ref = c.data(:,2); 

k = 1.4; %Specific heat ratio of air 
R = 8314; %Universal gas constant (J/kmolK)
M_air = 28.97; %Molar mass of air (kg/kmol)

M = v/(sqrt(k*R*T/M_air)); %Mach number
Re = rho * v * D_o / mu; %Reynolds number

Cd = interp1(M_ref, Cd_ref, M); %Drag coefficient dependent on drag coefficient. NOTE: With rocket aerodynamics simulation this reference data needs updating
D = 0.5 * Cd * rho * v^2 * pi/4 * D_o^2; %Drag (N)


