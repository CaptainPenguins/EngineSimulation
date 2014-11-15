function [mdot, P1] = orifice(PT, P2, D1, D2, N, L_p, Dh, rho, V, e) 

k = 1.31; %Specific heat ratio for N2O

A2 = N * pi/4 * D2^2; %Orifice Area

Cd = 0.29; %Discharge Coefficient. NOTE: This will change as the injector simulation becomes available
Beta = D2/D1;  %Ratio of hole diameter to pipe diameter
C = Cd / sqrt(1-Beta^4); %Flow Coefficient

[P1] = PipeLosses(PT, L_p, Dh, rho, V, e);  
r = P2/P1; 
Y = 1 - (1-r)/k * (0.41 + 0.35*Beta^4); 

mdot = C * A2 * Y * sqrt(2*rho*(P1-P2));

 






