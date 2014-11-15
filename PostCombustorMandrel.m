clear all; close all; 

r1 = (3.67+0.23)*25.4/2; %Outer radius (mm), 0.23 to account for the cut thickness
r2 = 1.5*25.4/2; %Mandrel radius (mm)

%%% Create circles %%%
theta = 0:0.01:2*pi;

x1 = r1 * cos(theta); x1 = [x1, x1(1)];
y1 = r1 * sin(theta); y1 = [y1, y1(1)]; 

x2 = r2 * cos(theta); x2 = [x2, x2(1)];
y2 = r2 * sin(theta); y2 = [y2, y2(1)]; 

%%% Create lead-in line %%%

x3 = x1(end): (x2(1) - x1(end))/10 :x2(1); 
y3 = zeros(1, length(x3));

%%% Stack contours %%%

x0 = [x1,x3,x2, fliplr(x3)]; 
y0 = [y1, y3, y2, fliplr(y3)]; 

%%% Reduce array size %%%

c = 1:1000; 
c0 = 1:1000/(length(x0)+1):1000; 

x = -interp1(c0, x0, c); 
y = interp1(c0, y0, c); 



