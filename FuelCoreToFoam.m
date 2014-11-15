clear all; close all; 

ri = 4.5 * 0.0254 / 2; %Radius of inner engine (m)
R_s = 0.44 * 0.0254 / 2; %Half-width of spoke (m)
N_s = 10; %Number of spokes
r_ic = 1.25 * 0.0254; %Radius of inner core (m)
r_s = 1.6 * 0.0254; %Radial distance to center of spoke
C = 0; %0 = spokes, 1 = circle
[x_core,y_core] = CoreShapeGenerator(ri, R_s, N_s, r_ic, r_s, C); 

%%% Create Mandrel Circle %%%

r_mc = 1.5 * 0.0254 / 2; %Radius of mandrel circle (m)
theta = 0:0.01:2*pi; 
x_mc = permute(r_mc * cos(theta), [2,1]);y_mc = permute(r_mc * sin(theta), [2,1]); 

%%% Create first lead-in line %%%

x_li1 = x_core(1): ((x_core(1)+0.5*0.0254) - (x_core(1)))/100 :x_core(1)+0.5*0.0254;
x_li1 = flipud(permute(x_li1, [2,1])); 
y_li1 = zeros(length(x_li1), 1); 

%%% Create second lead-in line %%%

x_li2 = x_core(end): (x_mc(1) - x_core(end))/100: x_mc(1); 
x_li2 = permute(x_li2, [2,1]); 
x_li2(end) = []; 
y_li2 = zeros(length(x_li2), 1); 

%%% Create lead-out line %%%

x_lo = x_mc(end): (x_li1(1) - x_mc(end))/100: x_li1(1);
x_lo = permute(x_lo, [2,1]); 
y_lo = zeros(length(x_lo), 1); 

%%% Combine arrays %%%

x0 = vertcat(x_li1, x_core, x_li2, x_mc, x_lo)*1000; 
y0 = vertcat(y_li1, y_core, y_li2, y_mc, y_lo)*1000; 

%%% Reduce array size %%%

c = 1:1500; 
c0 = 1:1500/(length(x0)+2):1500; 

x = -permute(interp1(c0, x0, c), [2,1]); x = x - min(x); 
y = permute(interp1(c0, y0, c), [2,1]); y = y - min(y); 

