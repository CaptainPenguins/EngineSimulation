function [mdotOx, mL, mv, m, T, P, burnout, Pi] = oxidizer (mv, T0, mL0, m0, Pcc, D_Ii, D_p, N_p, dt, V, L_p, Dh, mdotOx0, e)

burnout = 0; %Logical test to see if oxidizer is out. 0 = not out, 1 = out

Cv = 40; %Constant volume specific heat capacity (kJ/kmolK) NOTE: This value is relatively constant over our range of P and Ts

c = importdata('HeatOfVaporization.txt');
T_ref = c.data(:,1); Hvap_ref = c.data(:,2);

c = importdata('N2OThermophysics.txt');
T_rho = c.data(:,1); rhoV_ref = c.data(:,2); rhoL_ref = c.data(:,3); Pvap_ref = c.data(:,4);

Hvap = interp1(T_ref, Hvap_ref, T0)*1000; %Heat of vaporization (kJ/kmol)

Q = mv * Hvap; %Heat removed during vaporization (kJ)

T = -Q/(mL0 * Cv) + T0;    %New oxidizer tank temperature (K)
P = interp1(T_rho, Pvap_ref, T);  %Vapour pressure (Pa)

rhoL = interp1(T_rho, rhoL_ref, T); %Liquid density (kg/m^3)
rhoV = interp1(T_rho, rhoV_ref, T); %Vapour density (kg/m^3)

U = mdotOx0/(rhoL*pi/4 * Dh^2); %Flow velocity in plumbing

[mdotOx, Pi] = orifice(P, Pcc, D_Ii, D_p, N_p, L_p, Dh, rhoL, U, e);   %Calculate injector mass flow rate (kg/s)

mL_NE = mL0 - mdotOx*dt; %Mass of liquid in tank if N2O didn't react to expansion of nitrous vapour (kg)
m = m0 - mdotOx*dt; %New liquid + oxidizer mass (kg)

mL = (V - m/rhoV) / (1/rhoL - 1/rhoV);  %Actual mass of liquid left in tank (kg)

mv = mL_NE - mL; %Vapourized mass (kg)

if mv < 0
    mv = 0;
else if mL <= 0.1 %When mass of liquid is less than 50g, stop burn
        burnout = 1;
        mL = 0;
    end
end
