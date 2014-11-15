function [Pcc, Tcc, x, y, L, kP, MP, mdotF, mdotOx, POx, TOx, mOxL, mOx, mvOx, BO, Ac, OF, rdot, Pi] = Engine(D_th, D_i, D_Ii, D_p, N_p, x0, y0, L0, l_c, dt, r_mf, Pcc, POx, mdotOx, mv, TOx0, mOxL0, mOx0, V, Ac, rhoF, Tcc, L_p, Dh, e)

dPcc = 1e10; %Initialize, never used
dmdotOx = 1e10; %Initialize, never used

%% Constants %%
kF = 1.3; kOx = 1.4; %Specific heat ratios of fuel and oxidizer
R = 8314; %Universal gas constant (J/kmolK)
UR = 0.1; %Under-relaxation factor for chamber pressure
Tol = 1*6894.757; %Chamber pressure tolerance
OTol = 0.01; %Oxidizer mass flow rate tolerance

M_paraffin = 352.77; %molar mass of paraffin (kg/kmol)
M_Al = 26.98; %molar mass of aluminum (kg/kmol)
M_N2O = 44; %molar mass of nitrous (kg/kmol)
M_CO2 = 44.01; %Molar mass of CO2 (kg/kmol)
M_H2O = 18.01; %Molar mass of H2O (kg/kmol)
M_N2 = 28.02; %Molar mass of N2 (kg/kmol)
M_Al2O3 = 101.96; %Molar mass of Al2O3 (kg/kmol)

while (abs(dPcc) > Tol) || (abs(dmdotOx) > OTol)
        
    % Calculate oxidizer tank properties
    [mdotOxT, mOxL, mvOx, mOx, TOx, POx, BO, Pi] = oxidizer(mv, TOx0, mOxL0, mOx0, Pcc, D_Ii, D_p, N_p, dt, V, L_p, Dh, mdotOx, e);
    
    dmdotOx = mdotOxT - mdotOx; %Change in flow rate (kg/s)
    mdotOx = mdotOx + UR*dmdotOx; %New oxidizer mass flow rate (kg/s)
    Gox = mdotOx * 1000 / (Ac*100^2); %Mass flux (g/cm^2 s). 
   
    % Calculate regression rate
    rdot =  (0.155 *(1.1+1.25)/2 + 0.488)/2 * Gox ^ ((0.5+0.62)/2);  %Regression rate of R = 0.155Gox^0.5 * 17% (mm/s)
    
    % Calculate new fuel core geometry and arc length
    [x, y, L, Ac] = FuelCoreBurn(D_i/2, rdot/1000, dt, x0, y0);
    
    % Calculate combustion chamber kinematics
    mdotF = rhoF * l_c * 0.5 * (L+L0) * (rdot/1000);   %Fuel mass flow rate (kg/s)
    OF = mdotOx/mdotF;   %Oxidizer-Fuel Ratio
    kP = 1.15; %Value obtain from NASA CEA
 
    % Calculate Stoichiometrics
    A = [1 0 0 0 0 0 0 0;
        -M_Al*OF M_N2O 0 0 0 0 0 0;
        0 0 1 0 0 25 0 0;
        0 0 0 1 0 52 0 0;
        1 0 0 0 0 0 -1 -2;
        0 2 0 0 -2 0 0 0;
        0 1 -2 -1 0 0 0 -3;
        0 0 0 0 0 -r_mf*M_paraffin/M_Al 1 0];
    
    B = [r_mf*M_paraffin/M_Al; M_paraffin*OF; 25; 52; 0; 0; 0; 0];
    
    %% n(1) = Al
    %% n(2) = N2O
    %% n(3) = CO2
    %% n(4) = H2O
    %% n(5) = N2
    %% n(6) = C25H52
    %% n(7) = Al
    %% n(8) = Al2O3
    
    n = linsolve(A,B);
    nP =  n(3) + n(4) + n(5) + n(6) + n(7) + n(8); %Moles of products (mol)
    MP = 1/nP * (n(3)*M_CO2 + n(4)*M_H2O + n(5)*M_N2 + n(6)*M_paraffin + n(7)*M_Al + n(8)*M_Al2O3); %Molar mass of products (kg/kmol)
    
    % Calculate Combustion Chamber Temperature
    [Tcc] = ChamberTemperature(Pcc, n, Tcc);
   
    % Calculate Combustion Chamber Pressure
    cStar = sqrt(R*Tcc/(kP*MP) / ((2/(kP+1))^((kP+1)/(kP-1)))); %Characteristic exhaust velocity (m/s)
    Afc = (L+L0)/2 * l_c; %Burn surface area (m^2)
    AStar = pi/4 * D_th^2; %Nozzle throat area (m^2)
    dPcc = cStar / AStar * (Afc * rhoF * rdot/1000 + mdotOx) - Pcc; %Change in combustion chamber pressure (Pa)
    Pcc = Pcc + UR*dPcc; %Apply under-relaxation factor for combustion chamber pressure (Pa)
    
 
end

