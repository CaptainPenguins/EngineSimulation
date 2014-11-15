clear all; close all;

%% This is the main program for the 2014-15 Bia II engine simulation. A number of improvements have been made over last year's code, including:

% - Time stepped calculation of combusution chamber temperature
% - 2D regression-rate data for wagon wheel fuel core geometry
% - Engine geometry optimization algorithm (tentative as of Oct 24, 2013)
% - Model of oxidizer tank emptying
% - Comparison between 1st-order explicit time-marching method, and
% 2nd-order explicit and 1st-order implicit (tentative as of Oct 24, 2013)
% - Nozzle loss factors
% - Fluctuation in the paraffin-aluminum mixture composition to account for
% non-uniformities in the fuel core
% - Inclusion of launch angle into velocity calculations
% - Piping Loss Factors with coolant jacket losses
% - Increase in nozzle throat diameter due to erosion

fprintf('---------------  NEW ENGINE SIMULATION --------------\n');
fprintf('WARNING: THIS CODE REQUIRES ABOUT 5 HOURS TO COMPLETE\n'); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inputs: TMM, T_Ox, T_a, T_f, D_o, D_th, D_e, D_p, D_Ii, D_i
% Output: T_rho, rhoL_ref, rhoF, Ac0, CTr, m_F
% Functions called: CoreShapeGenerator.m, explicit.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

TMM = 1; %Time marching method, 1 = 1st-order explicit, 2 = 2nd-order explicit, 3 = implicit

%%% Ambient Conditions %%%
T_Ox = 293; %Oxidizer tank firing temperature (K)
T_a = 313.15; %Expected ambient temperature (K)
T_f = 273; %Ox tank filling temperature (K)

%%% Geometric Conditions to optimize %%%
D_o = 6 * 0.0254; %Rocket OD (m)
D_th = 1.1 * 0.0254; %Nozzle throat diameter (m)
D_e = 1.95 * 0.0254; %Nozzle exit diameter (m)
D_p = 7/64 * 0.0254; %Injector pore diameter (m)
D_Ii = 2.75 * 0.0254; %Injector inlet diameter (m)
D_i = 4.5 * 0.0254; %Engine inner diameter (m) 

%%% Fuel Core Shape Parameters %%%
W_s = 0.44 * 0.0254; %Width of spoke (m)
N_s = 10; %Number of spokes
R_ci = 1.25 * 0.0254; %Radius of inner core (m)
R_s = 1.6 * 0.0254; %Radius of spoke length (m)
C = 0; %0 = spokes, 1 = circle

N_p = 20; %Number of injector pores

l_c = 17.81 * 0.0254; %Length of fuel core (m)
l_d = 1.5 * 0.0254; %Length of divergent portion of nozzle (m)

r = 0.35; %Al-Paraffin mass fraction

%%% Oxidizer Tank Properties %%%
VOx = 5000* (1/100)^3; %Volume of flight tank (m^3)
FF = 0.95; %Fill factor (Vl/Vtot)

%%% Plumbing Loss Factors %%%
L_p = 3; %Length of plumbing (m)
Dh = 0.5 * 0.0254; %Hydraulic diameter of plumbing (m)
e = 2; %Number of 90deg bends in plumbing   

m0 = 37.4; %Take-off weight (kg)

[D, v, a, Pcc, Tcc, Cd, h, Th, M, I, OF] = EngineSimulation(TMM, T_Ox, T_a, D_o, D_th, D_e, D_p, D_Ii, D_i, N_p, l_c, r, W_s, N_s, R_ci, R_s, VOx, FF, l_d, T_f, L_p, Dh, e, C, m0);

fprintf('Altitude: %0.1f ft AGL\n', h(end)*3.28); 
fprintf('Max Mach: %0.2f\n', M); 
fprintf('Average O/F: %0.1f\n', OF); 
fprintf('Total Impulse: %0.1f Ns\n', I); 

%% Create Plots %%




