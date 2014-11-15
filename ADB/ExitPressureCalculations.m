function [Pn] = ExitPressureCalculations (K, At, Ax, Pc, Pa)

%Set up inital values and variables

P_exitGuess = Pa;
P_exitTemp = Pc*((K-1)/(K+1))*(1/(1-P_exitGuess/Pc)^((K-1)/K))*(At/(Ax*((K+1)/2)^(1/(K-1))))^(2*K);
tolerance = 0.01;

while (P_exitTemp - P_exitGuess)/P_exitGuess > tolerance
    
    P_exitGuess = P_exitTemp;
    P_exitTemp = Pc*((K-1)/(K+1))*(1/(1-P_exitGuess/Pc)^((K-1)/K))*(At/(Ax*((K+1)/2)^(1/(K-1))))^(2*K);
end

Pn = P_exitTemp;

%
%
% i = 2;
% Px = zeros(10000, 1);
% Px(1) = Pa;
%
%     Px(i) = Pc*((K-1)/(K+1))*(1/(1-Px(i-1)/Pc)^((K-1)/K))*(At/(Ax*((K+1)/2)^(1/(K-1))))^(2*K);
%     Pn = Px(i);
%     i = i+1;
%
% while(Px(i)-Px(i-1)>0.000001)
%
%     Px(i) = Pc*((K-1)/(K+1))*(1/(1-Px(i-1)/Pc)^((K-1)/K))*(At/(Ax*((K+1)/2)^(1/(K-1))))^(2*K);
%     Pn = Px(i);
%     i = i+1;
% end

