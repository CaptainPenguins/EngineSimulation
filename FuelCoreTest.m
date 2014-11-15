clear all; close all;

r_dot = 5/1000;
dt = 0.1;
ri = 2.253 * 0.0254; %Engine inner diameter (m)

% fBurn = '\\SRVA\Homes$\debiasia\Desktop\FuelCoreBurn\CorePics'; 

%%% CoreShapeGenerator Inputs (all in metric): 
% [x,y] = CoreShapeGenerator(<CC IR>, <Spoke Width>, <# of spokes>,<Core Inner Radius>, <Radius of spoke>);  
[x,y] = CoreShapeGenerator(4.5/2*0.0254, 0.0254*0.22, 10, 0.0254*1.25, 0.0254*1.6, 0); 
x0 = x; y0 = y;

phi = 0:0.01:2*pi; 
x_wall = ri*cos(phi); 
y_wall = ri*sin(phi); 

i = 1; 

%%% Initialize the Arc Length variable %%%
L(i) = 0; 
for j = 2:length(x0)
    L(i) = L(i) + sqrt((x0(j) - x0(j-1))^2 + (y0(j) - y0(j-1))^2); %Based only on initial profile
end

while L(i) > 0
    i = i+1; 
    [x, y, L(i)] = FuelCoreBurn(ri, r_dot, dt, x, y); 
    plot(x0/0.0254, y0/0.0254,'b', x/0.0254 ,y/0.0254, 'r', x_wall/0.0254, y_wall/0.0254, 'k'); axis equal; xlabel('x-position (in)'); ylabel('y-position (in)'); 
    xlim([-ri*1.25/0.0254 ri*1.25/0.0254]); ylim([-ri*1.25/0.0254 ri*1.25/0.0254]);
    h_legend = legend('Initial Profile', 'Burnt Profile', 'Combustion Chamber Wall'); 
    set(gcf,'units','pixel');set(gcf,'position',[0,0,960,960]);
    set(gcf,'papersize',[960,960]);
    Time = text(-3,2.5,['Time = ', num2str(i*dt),'s']); set(Time, 'FontSize',15); 
 pause(0.001); 
end

figure, plot((1:i)*dt, L*25.4); xlabel('Time (s)'); ylabel('Burn Arc Length (inch)'); 





