function [Pcc, OMFlow] = CombustionChamberProperties(CCPressure, OxPressure, At, K, r, CCTemp,densityOx,densityGrain,BArea,OrRad,OrNum,TotalOrArea,N,A,Beta,c)

%INITIALIZING INITIAL CONDITIONS

tolerance = 0.0001;

Pcc_Guess = CCPressure;
%ActualPcc = Pcc_Guess;

rr = Pcc_Guess/OxPressure;
Y = 1-((1-rr)/K)*(0.41 + 0.35*Beta); % Compressibility Factor

OMFlow_Guess = c*Y*TotalOrArea*sqrt(2*(OxPressure-Pcc_Guess)*densityOx);
OMFlux = OMFlow_Guess/((pi()*(OrRad)^2)*OrNum);
FuelRegressionRate = A*OMFlux^N;
Pcc_temp = (OMFlow_Guess + BArea*FuelRegressionRate*(densityGrain))/(At*sqrt(K/(r*CCTemp))*(2/(K+1))^((K+1)/(2*(K-1)))); %Calculates Combustion Chamber Pressure
ActualPcc = Pcc_Guess + 0.1*(Pcc_temp - Pcc_Guess);

%CALCULATING COMBUSTION CHAMBER CONDITIONS
while (ActualPcc - Pcc_Guess)/Pcc_Guess > tolerance
    Pcc_Guess = ActualPcc;
    
    rr = Pcc_Guess/OxPressure;
    %beta = (2*OrificeRad)/PipeDiameter; % Ratio of Orifice Area to Pipe Area
    Y = 1-((1-rr)/K)*(0.41 + 0.35*Beta); % Compressibility Factor
    
    OMFlow_Guess = c*Y*TotalOrArea*sqrt(2*(OxPressure-CCPressure)*densityOx);    
    OMFlux = OMFlow_Guess/((pi()*(OrRad)^2)*OrNum);

    FuelRegressionRate = A*OMFlux^N;

    Pcc_temp = (OMFlow_Guess + BArea*FuelRegressionRate*(densityGrain))/(At*sqrt(K/(r*CCTemp))*(2/(K+1))^((K+1)/(2*(K-1)))); %Calculates Combustion Chamber Pressure
    ActualPcc = Pcc_Guess + 0.1*(Pcc_temp - Pcc_Guess);
end

Pcc = ActualPcc;
OMFlow = OMFlow_Guess;