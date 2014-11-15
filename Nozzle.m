function [Th, Me, Pe, mdotT, Ve, Te] = Nozzle(P0, T0, k, D_th, D_e, M, Pa, mdotF, mdotOx, Mg, l_d, CTr)

R = 8314; %Universal gas constant (J/kmolK)
alpha = atan(((D_e - D_th)/2)/l_d); %Divergence angle (rad)
lambda = 0.5 * (1 + cos(alpha)); %Divergence correction factor

c = importdata('ChamberToThroat.txt'); 
CTr_ref = c.data(:,1); fTP_ref = c.data(:,2); 
fTP = interp1(CTr_ref, fTP_ref, CTr)/100; 

Pt = P0 / (1 + (k-1)/2)^(k/(k-1)) * fTP; 

% (1 + 0.5 * (k-1)^(k/(k+1)) * P0) * fTP; %Throat pressure (Pa)

tol = 1e-3; %Exit mach number tolerance

Me = (((D_e^2 * Mg / D_th^2)^((2*k-2)/(k+1)) * (1+(k-1)/2) - 1) /((k-1)/2))^0.5; %Exit mach number

while abs((Me - Mg)/Mg) > tol
    Mg = Me;
    Me = (((D_e^2 * Mg / D_th^2)^((2*k-2)/(k+1)) * (1+(k-1)/2) - 1) /((k-1)/2))^0.5;  %Exit mach number
end

Pe = Pt * ((1+ (k-1)/2 * Me^2)^(-k/(k-1))); 
Te = ((Pe/P0)^((k-1)/k) * T0)*0.995; %Exit Temperature with losses from chemical reactions in nozzle(K)
Ve = Me * sqrt(k*R/M*Te) * lambda * (1-0.75/100); %Exit velocity with divergence angle correction factor and wall friction factor (m/s)
mdotT = P0 * pi/4 * D_th^2 * sqrt(k/(R/M*T0) * (2/(k+1))^((k+1)/(k-1)));  %Choked mass flow rate

Th = ((mdotF + mdotOx) * Ve + (Pe - Pa) * pi/4 * D_e^2)*0.95;%Thrust with losses associated with solid particles in exhaust (N)







