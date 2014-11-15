function [P] = PipeLosses(Pt, L, Dh, rho, V, e)

lambda = 0.025; %Estimated based on overall Reynolds number and friction factor

epsilon1 = 1.5; %minor loss coefficient for 90deg elbow
epsilon2 = 0; %minor loss coefficient for tee
epsilon3 = 0.05; %minor loss coefficient for ball valve
epsilon4 = 0; %minor loss coefficient for coolant jacket

P = Pt - (e*epsilon1+epsilon2+epsilon3+epsilon4)*V^2*rho/2 - lambda * (L/Dh) * rho * V^2/2;
