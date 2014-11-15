function [Pcc, Tcc, x, y, L, kP, MP, mdotF, mdotOx, POx, TOx, mOxL, mOx, mvOx, BO, Ac, OF, rdot, Pi] = Engine(D_th, D_i, D_Ii, D_p, N_p, x0, y0, L0, l_c, dt, r_mf, Pcc, POx, mdotOx, mv, TOx0, mOxL0, mOx0, V, Ac, rhoF, Tcc, L_p, Dh, e)

dPcc = 1e10; %Initialize, never used
dmdotOx = 1e10; %Initialize, never used

%% Constants %%
kF = 1.3; kOx = 1.31; %Specific heat ratios of fuel and oxidizer
R = 8314; %Universal gas constant (J/kmolK)
UR = 0.1; %Under-relaxation factor for chamber pressure
Tol = 1*6894.757; %Chamber pressure tolerance
OTol = 0.01; %Oxidizer mass flow rate tolerance

while (abs(dPcc) > Tol) || (abs(dmdotOx) > OTol)  
    
    % Calculate oxidizer tank properties
    [mdotOxT, mOxL, mvOx, mOx, TOx, POx, BO, Pi] = oxidizer(mv, TOx0, mOxL0, mOx0, Pcc, D_Ii, D_p, N_p, dt, V, L_p, Dh, mdotOx, e);
    
    dmdotOx = mdotOxT - mdotOx; %Change in flow rate (kg/s)
    mdotOx = mdotOx + UR*dmdotOx; %New oxidizer mass flow rate (kg/s)
    Gox = mdotOx * 1000 / (Ac*100^2); %Mass flux (g/cm^2 s). 
   
    % Calculate regression rate
    rdot =  0.155*Gox^0.5;  % OLD REGRESSION RATE FOR PARAFFIN-AL: (0.155 *(1.1+1.25)/2 + 0.488)/2 * Gox ^ ((0.5+0.62)/2) (mm/s)
    
%------------------------------------------------------------------------------------------------------------------------------------------------------------------------%    
    
    % Calculate new fuel core geometry and arc length
    [x, y, L, Ac] = FuelCoreBurn(D_i/2, rdot/1000, dt, x0, y0);
    
    % Calculate combustion chamber kinematics
    mdotF = rhoF * l_c * 0.5 * (L+L0) * (rdot/1000);   %Fuel mass flow rate (kg/s)
    OF = mdotOx/mdotF;   %Oxidizer-Fuel Ratio
     
    % Calculate Tcc(Pcc, OF), MP(Pcc, OF), kP(Pcc, OF) from degree = 3 polynomial regression   
    Tcc = 8.670444e+02 + 4.056212e-07*Pcc^3 + -4.842914e-05*Pcc^2*OF+2.116436e-02*Pcc*OF^2+-1.380426e+01*OF^3+-5.028846e-04*Pcc^2+-9.840913e-02*Pcc*OF+1.212366e+02*OF^2+4.406476e-01*Pcc+1.441427e+02*OF;
    MP = 1.037491e+01 + 1.166876e-09*Pcc^3 + 5.240254e-08*Pcc^2*OF+9.831748e-05*Pcc*OF^2+-2.313642e-02*OF^3+-2.209588e-06*Pcc^2+-8.768980e-04*Pcc*OF+2.069864e-02*OF^2+2.965474e-03*Pcc+3.168387e+00*OF;
    kP = 1.045825e+00 + 4.642170e-11*Pcc^3 + -1.295318e-08*Pcc^2*OF+-7.977625e-07*Pcc*OF^2+2.873440e-03*OF^3+-2.893506e-08*Pcc^2+2.569395e-05*Pcc*OF+-4.420558e-02*OF^2+-3.373320e-05*Pcc+1.797579e-01*OF;
    
    % Calculate Combustion Chamber Pressure
    cStar = sqrt(R*Tcc/(kP*MP) / ((2/(kP+1))^((kP+1)/(kP-1)))); %Characteristic exhaust velocity (m/s)
    Afc = (L+L0)/2 * l_c; %Area fuel core (m^2)
    AStar = pi/4 * D_th^2; %Nozzle throat area (m^2)
    
    dPcc = cStar / AStar * (Afc * rhoF * rdot/1000 + mdotOx) - Pcc; %Change in combustion chamber pressure (Pa) [SUSPECT! rhoF]
    Pcc = Pcc + UR*dPcc; %Apply under-relaxation factor for combustion chamber pressure (Pa)
    
 
end










