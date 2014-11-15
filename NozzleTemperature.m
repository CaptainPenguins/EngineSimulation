function [T] = NozzleTemperature(dx, M, h, Tcc, T0) 

K = importdata('GraphiteConductivity.txt'); 

T_ref = K.data(:,1); k_ref = K.data(:,2); 
k = interp1(T_ref, k_ref, T0); 

rho = 2150; %Density (kg/m^3)
Cp = 710; %Specific heat capacity (J/kgK)
nu = mean(k/(rho * Cp)); %thermal diffusivity (m^2/s)

A = 1/dx^2 * ((1*diag(ones(M-1, 1), -1) + diag(ones(M-1, 1), 1)) - 2*diag(ones(M,1), 0)); %Center-difference scheme
A(end, end) = -1/dx^2; %Apply the perfectly insulated BC
bc = zeros(M, 1);
bc(1,1) = Tcc/(dx^2);

T = (eye(M) + h * A * nu)*T0 + h*nu.*bc; %Explicit Euler time-marching method

%%% NOTE: Explicit time-marching is pretty unstable because of high
%%% gradients. Future codes will likely need an implicit method. 





