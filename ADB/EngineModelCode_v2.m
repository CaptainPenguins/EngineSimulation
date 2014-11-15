clear all; close all; tic; clc;

fprintf('\n-----------------------------\n NEW SIMULATION TEST\n-----------------------------\n'); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%SETING UP ARRAYS%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dt = 0.01; %In seconds, a small time interval
%burnTime = 4.5; % Sets the burn time of the engine

%{
TotalMass = zeros(burnTime/dt, 1); 
velocity = zeros(burnTime/dt, 1);  
displacement = zeros(burnTime/dt, 1); 
Temp = zeros(burnTime/dt, 1); 
AmbiantPressure = zeros(burnTime/dt, 1);
rho = zeros(burnTime/dt, 1); 
drag = zeros(burnTime/dt, 1); 
acceleration = zeros(burnTime/dt, 1); 
Time = dt:dt:burnTime;

CombustionChamberPressure = zeros(burnTime/dt, 1);
CombustionChamberTemperature = zeros(burnTime/dt, 1);
ActualCombustionChamberPressure = zeros(burnTime/dt, 1);
OxTankPressure = zeros(burnTime/dt, 1);
OxMassFlux = zeros(burnTime/dt, 1);
OxMassFlowRate = zeros(burnTime/dt, 1); 
BurnArea = zeros(burnTime/dt, 1);
FuelRegressionRate = zeros(burnTime/dt, 1);
ExitPressure = zeros(burnTime/dt, 1);
ExaustVelocity = zeros(burnTime/dt, 1);
Thrust = zeros(burnTime/dt, 1);

%Arrays to carry number of mols/molecule
AmolN2O = zeros(burnTime/dt, 1);
BmolAl = zeros(burnTime/dt, 1);
CmolCO2 = zeros(burnTime/dt, 1);
DmolH2O = zeros(burnTime/dt, 1);
EmolAl2O3 = zeros(burnTime/dt, 1);
FmolC25H52 = zeros(burnTime/dt, 1);
GmolAl = zeros(burnTime/dt, 1);
HmolN2 = zeros(burnTime/dt, 1);
%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%INICALIZE INPUT VARIABLES%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
i=1;
g = 9.81; %Gravitation Constant
Cd = 0.75;
TotalPropellantMass = 5.63; %The total mass of propellant in kg
diameter = 6 * 0.0254; %Rocket diameter in m
Angle = 85; %Launch angle
%avgMassFlow = TotalPropellantMass/burnTime; %The average mass flow rate
%Impulse = 12150; % Set total impulse of engine
%avgThrust = 2700;%Impulse/burnTime; % Calc avg thrust of engine

%rhoProducts = ; %[]
C=0.8;
rhoOx = 750; %[kg/m^3] Liquid phase
rhoGrain = 1250; %[kg/m^3]
%BurnArea(i) = 52.25 * 0.0254; %[m^2]
OrificeRad = 0.0625 * 0.0254; %[m]
OrificeNum = 15;
PipeDiameter = 0.5 * 0.0254; %[m]
TotalOrfaceArea = pi()*(OrificeRad^2)*OrificeNum; %[m^2]
n = 0.5;
%a = 0.1705; % 10% increase in performance
a = 0.182125; % 17.5% increase in performance
%a = 0.19375; % 0.18375% increase in performance
beta = (2*OrificeRad)/PipeDiameter; % Ratio of Orifice Area to Pipe Area

%Update in Pcc function!
ThroatArea = 0.66476 * 0.0254; %[m^2]
R = 8.314; %[J/molK] %0.08205746; %[L?atm?/K?mol]
To = 3500; %[K] - Tempurature in combustion chamber
NozzelExitArea = 2.1124 * 0.0254; %[m^2]
k = 1.35; %[]

r = 0.4/0.6; %Mass fraction of paraffin-aluminum
OtoF = 1; %Oxidizer-to-fuel ratio

r_dot = 5/1000;
dt = 0.1;
ri = 2.5235/25.4; %Engine IR
coreLength = 0.6096; %[m]
%%% CoreShapeGenerator Inputs (all in metric): 
% [x,y] = CoreShapeGenerator(<CC IR>, <Spoke Width>, <# of spokes>,<Core Inner Radius>, <Radius of spoke>);  
[x,y] = CoreShapeGenerator(ri, 0.5/25.4, 5, 1.15/25.4, 1.5/25.4); 
x0 = x; y0 = y;

phi = 0:0.01:2*pi; 
x_wall = ri*cos(phi); 
y_wall = ri*sin(phi); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%START UP CONDITIONS%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
displacement(1) = 1243; %1243 metres - Green River elavation.
TotalMass(1) = 30.4; %The total inital mass of the rocket in kg
[Temp(1), AmbiantPressure(1), rho(1)] = atmosphericConditions(displacement(1)); %Executes the atmosphericProperties function

%%% Initialize the Arc Length variable %%%
L(i) = 0; 
for j = 2:length(x0)
    L(i) = L(i) + sqrt((x0(j) - x0(j-1))^2 + (y0(j) - y0(j-1))^2); %Based only on initial profile
end

OxTankPressure(i) = 700*6894.75728; %Sets inital tank pressue to 700psi then converts to Pa.
CombustionChamberPressure(i) = 200*6894.75728 ; %Sets inital Chamber pressure to 200psi then converts to Pa.
[CombustionChamberTemperature(i)] = ChamberTemperature(r, OtoF, CombustionChamberPressure(i));
%ActualCombustionChamberPressure(1) = CombustionChamberPressure(1); %Sets correction array to 0 initally
%[OxMassFlowRate(1), OxMassFlux(1), BurnArea(1), FuelRegressionRate(1)] = CombustionChamberProperties(CombustionChamberPressure(1), OxTankPressure(1));%Executes the CombustionChamberProperties function

rr = CombustionChamberPressure(i)/OxTankPressure(i);
Y = 1-((1-rr)/k)*(0.41 + 0.35*beta); % Compressibility Factor
OxMassFlowRate(i) = C*Y*TotalOrfaceArea*sqrt(2*(OxTankPressure(i)-CombustionChamberPressure(i))*rhoOx);    
OxMassFlux(i) = OxMassFlowRate(i)/((pi()*(OrificeRad)^2)*OrificeNum);
FuelRegressionRate(i) = a*OxMassFlux(i)^n;

[ExitPressure(i)] = ExitPressureCalculations(k, ThroatArea, NozzelExitArea, CombustionChamberPressure(i), AmbiantPressure(1)); %Calculates the inital exit pressure with respect to the combustion chamber
ExaustVelocity(i) = sqrt((2*k)/(k-1)*R*To*(1-(ExitPressure(i)/CombustionChamberPressure(i))^((k-1)/k))); %Calculate the inital exaust velocity

drag(i) = 0.5 * rho(i) * Cd * velocity(i)^2 * pi/4 * diameter^2;
Thrust(i)  = (FuelRegressionRate(i)+OxMassFlowRate(i))*ExaustVelocity(i) + (ExitPressure(i) - AmbiantPressure(i))*NozzelExitArea;
acceleration(i) = (Thrust(i) - TotalMass(i) * g - drag(i,1))/TotalMass(i);

[AmolN2O(i), BmolAl(i), CmolCO2(i), DmolH2O(i), EmolAl2O3(i), FmolC25H52(i), GmolAl(i), HmolN2(i)] = ChemicalEquationBalance(OxMassFlowRate(i));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%BurnTime Flight%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
while L(i) > 0 %Executes until the engine stops burning
    i = i+1;
    TotalMass(i) = TotalMass(i-1) - dt * (FuelRegressionRate(i-1)*burnArea(i));
    velocity (i) = velocity(i-1)  + acceleration(i-1) * dt; 
    displacement(i) = displacement(i-1) + (velocity(i)+velocity(i-1))/2*dt;  
    
    [x, y, L(i)] = FuelCoreBurn(ri, r_dot, dt, x, y); 
    plot(x0*25.4, y0*25.4,'b', x*25.4 ,y*25.4, 'r', x_wall*25.4, y_wall*25.4, 'k'); axis equal; xlabel('x-position (in)'); ylabel('y-position (in)'); 
    xlim([-ri*1.25*25.4 ri*1.25*25.4]); ylim([-ri*1.25*25.4 ri*1.25*25.4]);
    h_legend = legend('Initial Profile', 'Burnt Profile', 'Combustion Chamber Wall'); 
    set(gcf,'units','pixel');set(gcf,'position',[0,0,960,960]);
    set(gcf,'papersize',[960,960]);
    Time = text(-3,2.5,['Time = ', num2str(i*dt),'s']); set(Time, 'FontSize',15); 
    pause(0.001);
    burnArea(i) = L(i) * coreLength;

    [Temp(i), AmbiantPressure(i), rho(i)] = atmosphericConditions(displacement(i)); 
    drag(i) = 0.5 * rho(i) * Cd(1) * velocity(i)^2 * pi/4 * diameter^2;
    
    [CombustionChamberPressure(i), OxMassFlowRate(i)] = CombustionChamberProperties(CombustionChamberPressure(i-1), OxTankPressure(i-1), ThroatArea, k, R, CombustionChamberTemperature(i-1),rhoOx,rhoGrain,burnArea(i),OrificeRad,OrificeNum,TotalOrfaceArea,n,a,beta,C);
    OxTankPressure(i) = 700*6894.75728;
    [AmolN2O(i), BmolAl(i), CmolCO2(i), DmolH2O(i), EmolAl2O3(i), FmolC25H52(i), GmolAl(i), HmolN2(i)] = ChemicalEquationBalance(OxMassFlowRate(i));
    [CombustionChamberTemperature(i)] = ChamberTemperature(r, OtoF, CombustionChamberPressure(i));
    
    [ExitPressure(i)] = ExitPressureCalculations(k, ThroatArea, NozzelExitArea, ActualCombustionChamberPressure(i), AmbiantPressure(i));
    ExaustVelocity(i) = sqrt((2*k)/(k-1)*R*To*(1-(ExitPressure(i)/CombustionChamberPressure(i))^((k-1)/k)));
 
    Thrust(i)  = (FuelRegressionRate(i)+OxMassFlowRate(i))*ExaustVelocity(i) + (ExitPressure(i) - AmbiantPressure(i))*NozzelExitArea;
    acceleration(i) = (Thrust(i) - TotalMass(i) * g - drag(i,1))/TotalMass(i); 
    
end

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %Unaided Assent Flight%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% while displacement(i) > displacement(i-1)
%     i = i+1;
%     
%     TotalMass(i) = TotalMass(i-1);
%     velocity (i) = velocity(i-1)  + acceleration(i-1) * dt; 
%     displacement(i) = displacement(i-1) + (velocity(i)+velocity(i-1))/2*dt;    
% 
%     [Temp(i), AmbiantPressure(i), rho(i)] = atmosphericConditions(displacement(i)); 
%     drag(i) = 0.5 * rho(i) * Cd(1) * velocity(i)^2 * pi/4 * diameter^2;
%     
%     acceleration(i) = (-TotalMass(i) * g - drag(i))/TotalMass(i);
%     
%     Time(i) = Time(i-1) + dt; %Update the time array
% 
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Freefall Flight%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
while displacement(i)>1243
    i = i+1;
    
    TotalMass(i) = TotalMass(i-1);
    velocity (i) = velocity(i-1)  + acceleration(i) * dt; 
    displacement(i) = displacement(i-1) + (velocity(i)+velocity(i-1))/2*dt;    

    [Temp(i), Pressure(i), rho(i)] = atmosphericConditions(displacement(i)); 
    drag(i) = 0.5 * rho(1) * Cd(1) * velocity(1)^2 * pi/4 * diameter^2;
    
    acceleration(i) = (mass(i) * g - drag(i,1))/mass(i);
    
    Time(i) = Time(i-1) + dt; %Update the time array
end
%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CREATE PLOTS %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[AX,H1,H2] = plotyy(Time, (displacement - displacement(1)) * 3.28, Time, velocity);
legend('Rocket Height', 'Rocket Speed'); 
set(get(AX(1),'Ylabel'),'String','Height (ft)') 
set(get(AX(2),'Ylabel'),'String','Velocity (m/s)')
xlabel('Time (sec)'); 

figure, plot((1:i)*dt, L*25.4); xlabel('Time (s)'); ylabel('Burn Arc Length (inch)'); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%RESULT STATEMENTS%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('\nSIMULATION FINISHED, Total Computation time = %0.1f seconds\n\n', toc); 
fprintf('\nMax Altitude Achieved = %0.0f ft AGL\n', (displacement(i) - displacement(1) ) * 3.28); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
