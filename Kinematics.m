function [a, v, h] = Kinematics(m, m0, T, D, v0, h0, dt)

phi = 85; %Launch angle

%% Constants

g = 9.81; %Gravitational constant
a = (T - D - ((m+m0)/2)*g)/((m+m0)/2); %Average acceleration (m/s^2)
v = a*dt + v0; %Final velocity (m/s)
h = h0 + (v0*dt + a*dt^2/2)*sin(phi*pi/180); %Final height (m)

