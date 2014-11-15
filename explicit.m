function [D, v, a, Pcc, Tcc, Cd, Dh, Th, MM, I, MOF] = explicit(ri, D_th0, D_e, T_a0, D_o, D_Ii, D_p, D_i, N_p, l_c, r0, m0, TOx0, VOx, FF, x, y, Ac0, rhoF, mF0, l_d, CTr, T_f, L_p, Dh, e)
%%list out variables in xlsx file


%% Set initial conditions %%
dt = 0.001; %I have found that if dt>~0.001 the code will crash. Total simulation time is around 5 hours
i = 1;
h(i) = 1318.3; %Initial launch height (m)
v(i) = 0;
PccG = 250 * 6894.76; %Initial combustion chamber pressure guess (Pa)
POxG = 700 * 6894.76; %Initial Oxidizer Tank Guess (Pa)
mdotOxG = 0.5; %Initial oxidizer mass flow rate guess (kg/s)
MeG = 2.5; %Exit mach number guess
t(i) = 0; %Time counter (s)
BO = 0; %Logical test of burnout
I = 0; %Total impulse
m(i) = m0;
mF(i) = mF0;
fD_th = 6/10; % Percent increase in throat diameter per second (%/s)

%%% Calculate Initial Oxidizer Mass %%%
c = importdata('N2OThermophysics.txt');
T_rho = c.data(:,1); rhoL_ref = c.data(:,3); rhoV_ref = c.data(:,2);
rhoL0 = interp1(T_rho, rhoL_ref, TOx0);
rhoV0 = interp1(T_rho, rhoV_ref, TOx0);
rhoLf = interp1(T_rho, rhoL_ref, T_f); %Filling density
rhoVf = interp1(T_rho, rhoL_ref, T_f);
mOx0 = VOx*FF*rhoLf + (1-FF)*VOx*rhoVf; %Initial total oxidizer mass (kg)
mOxL0 = VOx*FF*rhoLf; %Initial oxidizer liquid mass (kg)
mv0 = 0.01; %Initial guess for oxidizer vapourized mass (kg)
TccG = 3500; %Guess for adiabatic flame temperature

%%% Create geometry for nozzle temperature model %%%
t_N = 0.5 * 0.0254; %Thinnest graphite portion (m)
N = 5; %Number of nodes
T0 = T_a0 * ones(N,1); %Initial temperature profile

%%% Create a wall geometry for plotting %%%
x0 = x; y0 = y;
phi = 0:0.01:2*pi;
x_wall = ri*cos(phi);
y_wall = ri*sin(phi);

%%% Initialize the Arc Length variable %%%
L(i) = 0;
for j = 2:length(x0)
    L(i) = L(i) + sqrt((x0(j) - x0(j-1))^2 + (y0(j) - y0(j-1))^2); %Based only on initial profile
end

%%% Execute Burn Loop %%%

while L(i) > 0 && BO == 0
    
    r = r0 + (rand(1,1) - 0.5)*0.05; %Add a 5% fluctuating component to the paraffin-aluminum composition
    D_th(i) = D_th0*(1+fD_th/100*i*dt); %Accounts for erosion of nozzle throat
    
    % Calculate Ambient Conditions
    [T_a(i) P_a(i) rho_a(i), mu_a(i)] = atmosphericProperties(h(i), T_a0, h(1));
    
    % Calculate Aerodynamic Properties
    [Cd(i) D(i) M(i) Re(i)] = Aerodynamics(v(i), T_a(i), mu_a(i), rho_a(i), D_o);
    
    if i == 1
        % Calculate Engine Properties
        [Pcc(i), Tcc(i), x, y, L(i+1), kP, MP, mdotF(i), mdotOx(i), POx(i), TOx, mOxL, mOx(i), mv, BO, Ac, OF(i), rdot(i), Pi] = Engine(D_th(i), D_i, D_Ii, D_p, N_p, x, y, L(i), l_c, dt, r, PccG, POxG, mdotOxG, mv0, TOx0, mOxL0, mOx0, VOx, Ac0, rhoF, TccG, L_p, Dh, e);
        % Calculate Nozzle Temperature
        [Tn] = NozzleTemperature(t_N/N, N, dt, Tcc(i), T0);
        
        % Calculate Nozzle Conditions
        [Th(i), Me(i), Pe(i), mdotT(i), Ve(i), Te(i)] = Nozzle(Pcc(i), Tcc(i), kP, D_th(i), D_e, MP, P_a(i), mdotF(i), mdotOx(i), MeG, l_d, CTr);
        I = Th(i)*dt + I; %Calculate total impulse
    else
        % Calculate Engine Properties, includes increase of throat diameter
        % as a result of erosion
        [Pcc(i), Tcc(i), x, y, L(i+1), kP, MP, mdotF(i), mdotOx(i), POx(i), TOx(i), mOxL, mOx(i), mv, BO, Ac, OF(i), rdot(i), Pi] = Engine(D_th(i), D_i, D_Ii, D_p, N_p, x, y, L(i), l_c, dt, r, Pcc(i-1), POx(i-1), mdotOx(i-1), mv, TOx(i-1), mOxL, mOx(i-1), VOx, Ac, rhoF, Tcc(i-1), L_p, Dh, e);
        
        % Calculate Nozzle Temperature
        [Tn] = NozzleTemperature(t_N/N, N, dt, Tcc(i), Tn);
  
            
        % Calculate Nozzle Conditions, includes increase of throat diameter
        % as a result of erosion
        [Th(i), Me(i), Pe(i), mdotT(i), Ve(i), Te(i)] = Nozzle(Pcc(i), Tcc(i), kP, D_th(i), D_e, MP, P_a(i), mdotF(i), mdotOx(i), Me(i-1), l_d, CTr);
        I = Th(i)*dt + I; %Calculate total impulse
    end
    
    % Calculate new masses
    mF(i+1) = mF(i) - mdotF(i)*dt; %(kg)
    m(i+1) = m(i) - (mdotF(i) + mdotOx(i))*dt; %(kg)
    
    % Calculate kinematics
    [a(i), v(i+1), h(i+1)] = Kinematics(m(i+1), m(i), Th(i), D(i), v(i), h(i), dt);
    t(i+1) = t(i) + dt;
    
    i = i+1;   %Step counter
    plot(x0/0.0254, y0/0.0254,'b', x/0.0254 ,y/0.0254, 'r', x_wall/0.0254, y_wall/0.0254, 'k'); axis equal; xlabel('x-position (in)'); ylabel('y-position (in)');
    xlim([-ri*1.25/0.0254 ri*1.25/0.0254]); ylim([-ri*1.25/0.0254 ri*1.25/0.0254]);
    legend('Initial Profile', 'Burnt Profile', 'Combustion Chamber Wall');
    Time = text(-2.75,2.5,['Time = ', num2str(t(i)),'s']); set(Time, 'FontSize',10);
    set(gcf,'units','pixel');set(gcf,'position',[0,0,600,600]);
    set(gcf,'papersize',[600,600]);
    pause(0.001);
    
    
    
    
end

D_th(i) = D_th(i-1);
Te(i) = Te(i-1);
Me(i) = 0;
Ve(i) = 0;
D(i) = D(i-1);
Th(i) = 0;
Pe(i) = Pe(i-1);
P_a(i) = P_a(i-1);
POx(i) = P_a(i-1); Pcc(i) = P_a(i-1);
TOx(i) = TOx(i-1); Tcc(i) = Tcc(i-1);
mOx(i) = mOx(i-1);
OF(i) = OF(i-1);
rdot(i) = 0;

% Plot Engine and Oxidizer Tank Properties %
figure;
subplot(2,2,1);
plot(t, POx/6894, t, Pcc/6894, t, Pe/6894, t, P_a/6894);
xlabel('Time (s)'); ylabel('Pressure (psi)');
legend('Oxidizer Tank', 'Combustion Chamber', 'Nozzle Exit', 'Ambient');

subplot(2,2,2);
[AX, H1, H2] = plotyy(t, TOx - 273, t, Tcc - 273);
set(get(AX(1), 'Ylabel'), 'String', 'Oxidizer Tank Temperature (C)', 'FontSize', 10);
set(get(AX(2), 'Ylabel'), 'String', 'Combustion Chamber Temperature (C)', 'FontSize', 10);
xlabel('Time (s)');

subplot(2,2,3);
plot(t, m, 'r', t, mF, 'b', t, mOx, 'g'); xlabel('Time (s)'); ylabel('Mass (kg)');
legend('Total Mass', 'Fuel Mass', 'Oxidizer Mass');

subplot(2,2,4);
[AX, H1, H2] = plotyy(t, OF, t, rdot); xlabel('Time (s)'); ylabel('Oxidizer-Fuel Ratio');
set(get(AX(1), 'Ylabel'), 'String', 'Oxidizer-Fuel Ratio');
set(get(AX(2), 'Ylabel'), 'String', 'Regression Rate (mm/s)');
xlabel('Time (s)');

set(gcf,'units','pixel');set(gcf,'position',[0,0,600,600]);
set(gcf,'papersize',[800,800]);

figure;
subplot(2,2,1);
plot(t, Th, t, D); xlabel('Time (s)'); ylabel('Force (N)'); legend('Thrust', 'Drag');

subplot(2,2,2);
[AX, H1, H2] = plotyy(t, Ve, t, Me);
set(get(AX(1), 'Ylabel'), 'String', 'Exhaust Velocity (m/s)');
set(get(AX(2), 'Ylabel'), 'String', 'Exhaust Mach Number');
xlabel('Time (s)');

subplot(2,2,3);
[AX, H1, H2] = plotyy(t, Te - 273, t, D_th*39.37);
set(get(AX(1), 'Ylabel'), 'String', 'Exhaust Temperature (C)');
set(get(AX(2), 'Ylabel'), 'String', 'Nozzle Throat Diameter (inch)');
xlabel('Time (s)');

subplot(2,2,4);
plot(t(1:i-1), mdotT, t(1:i-1), mdotF, t(1:i-1), mdotOx); xlabel('Time (s)'); ylabel('Mass flow rate (kg/s)');
legend('Choked Flow Rate', 'Fuel Mass Flow Rate', 'Oxidizer Mass Flow Rate');

set(gcf,'units','pixel');set(gcf,'position',[0,0,600,600]);
set(gcf,'papersize',[800,800]);

% Export data to Excel spreadsheet
% T = table(t, Pcc/6894, Tcc - 273, OF, rdot); % Time (s), Combustion Chamber Pressure (psia), Combustion Chamber Temprature (C), Oxidiser-Fuel Ratio, Regression Rate (mm/s)
% filename = 'testrundata.xlsx';
% writeTable(T, filename);

% figure; 
% plot(0+(t_N/N)/2:t_N/N:t_N - (t_N/N)/2, Tn-273); xlabel('X-position (m)'); ylabel('Temperature (C)'); 

%%% Execute Free Flight Loop %%%
while v(i) > 0
    
    Th(i) = 0; %Thrust = 0
    
    % Calculate Ambient Conditions
    [T_a(i) P_a(i) rho_a(i), mu_a(i)] = atmosphericProperties(h(i), T_a0, h(1));
    
    % Calculate Aerodynamic Properties
    [Cd(i) D(i) M(i) Re(i)] = Aerodynamics(v(i), T_a(i), mu_a(i), rho_a(i), D_o);
    
    % Update the mass array
    m(i+1) = m(i);
    
    % Calculate kinematics
    [a(i), v(i+1), h(i+1)] = Kinematics(m(i+1), m(i), Th(i), D(i), v(i), h(i), dt);
    t(i+1) = t(i) + dt;
    
    i = i+1;   %Step counter
    
end

% Create variables to return for print statements
Dh = h(end) - h(1); %Total change in altitude
MM = max(M); %Maximum mach number
MOF = mean(OF); %Mean O/F

% Plot Trajectory %
figure;
[ax, h1, h2] = plotyy(t,v, t, (h - h(1))*3.28);
set(get(ax(1), 'Ylabel'), 'String', 'Velocity (m/s)', 'FontSize', 10);
set(ax(1), 'YLim', [0 400]);
set(ax(1), 'YTick', [0:50:400]);
set(get(ax(2), 'Ylabel'), 'String', 'Height (ft)', 'FontSize', 10);
set(ax(2), 'YLim', [0 15000]);
set(ax(2), 'YTick', [0:2500:15000]);
set(gcf,'units','pixel');set(gcf,'position',[0,0,600,600]);
set(gcf,'papersize',[600,600]);














