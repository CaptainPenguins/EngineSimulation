function [D, v, a, Pcc, Tcc, Cd, h, Th, M, I, OF] = EngineSimulation(TMM, T_Ox, T_a, D_o, D_th, D_e, D_p, D_Ii, D_i, N_p, l_c, r, W_s, N_s, R_ci, R_s, VOx, FF, l_d, T_f, L_p, Dh, e, C, m0)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inputs: N2O thermophysical properties, rhoP, rhoAl
% Output: T_rho, rhoL_ref, rhoF, Ac0, CTr, m_F
% Functions called: CoreShapeGenerator.m, explicit.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% IMPORT VARIABLES %%

c = importdata('N2OThermophysics.txt');
T_rho = c.data(:,1); rhoL_ref = c.data(:,3);
rhoP = 900; 
rhoAl = 2700;
rhoF = 1/(r/rhoAl + (1-r)/rhoP); %Fuel density

%%% Create Initial Core Shape Geometry %%%
[x,y] = CoreShapeGenerator(D_i/2, W_s/2, N_s, R_ci, R_s, C);
Ac0 = polyarea(x,y); %Open area (m^2)

fprintf('Minimum Combustor Cross-Sectional Area: %0.5f\n', Ac0); 

CTr = Ac0 / (pi/4 * D_th^2); %Chamber-to-Throat area ratio

%% INITIAL CONDITIONS %%
m_F = (pi/4*D_i^2 - Ac0)*l_c*rhoF; %Mass of fuel

%% SELECT TIME MARCHING METHOD %%

switch TMM %[SUSPECT! Implement other time-marching methods; investigate scheduled relaxation Jacobi method?]
    case 1 
        [D, v, a, Pcc, Tcc, Cd, h, Th, M, I, OF] = explicit(D_i/2, D_th, D_e, T_a, D_o, D_Ii, D_p, D_i, N_p, l_c, r, m0, T_Ox, VOx, FF, x, y, Ac0, rhoF, m_F, l_d, CTr, T_f, L_p, Dh, e);
end