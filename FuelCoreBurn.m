function [x_new, y_new, L, A] = FuelCoreBurn(ri, r_dot, dt, x0, y0)

%% Author: J. Osborne, October 2013

% This code will take the initial (x,y) coordinates of a generic fuel core
% and, using regression rate data and time-stepping information, propagate
% the fuel core in a direction normal to the surface of the fuel core. The
% burn will continue until the chamber wall has been reached. The
% parameters for this code are:

%ri = Engine inside radius (m)
%r_dot = Regression rate (m/s)
%dt = Time step increment (s)
%x0 = original x-coordinates (m)
%y0 = original y-coordinates (m)
%x_new = final x-coordinates (m)
%y_new = final y-coordinates (m)
%L = Arc length of fuel core remaining (m)
%A = Open area of combustion chamber (m^2)

r0 = (x0.^2 + y0.^2).^0.5; 
theta0 = mod(atan2(y0, x0), 2*pi); 

theta = 0:0.01:2*pi; 
r = zeros(length(theta), 1); 
x = zeros(size(r)); y = zeros(size(r)); 

for i = 1:length(theta)
    r(i) = interp1(theta0, r0, theta(i), 'cubic'); 
    x(i) = r(i)*cos(theta(i)); 
    y(i) = r(i)*sin(theta(i)); 
end

%%%%% EXPAND CURVE TO NEXT TIME STEP %%%%%

dy = -gradient(x); dx = gradient(y);

delta_y = zeros(length(y), 2);
delta_x = zeros(length(x), 2);

x_new = zeros(length(x), 1);
y_new = zeros(length(x), 1);
theta_new = zeros(size(x_new)); 

for i = 1:length(x)
    
    delta_y(i,:) = [sqrt((r_dot*dt)^2 / ((dx(i)/dy(i))^2 + 1)), -sqrt((r_dot*dt)^2 / ((dx(i)/dy(i))^2 + 1))];
    delta_x(i,:) = [sqrt((r_dot*dt)^2 - delta_y(i,1)^2), -sqrt((r_dot*dt)^2 - delta_y(i,1)^2)];
    
    if sign(dx(i)) < 0
        
        x_new(i) = x(i) + delta_x(i,2);
        
    else if sign(dx(i)) > 0
            
            x_new(i) = x(i) + delta_x(i,1);
        else
            x_new(i) = x(i);
        end
    end
    
    if sign(dy(i)) < 0
        
        y_new(i) = y(i) + delta_y(i,2);
        
    else if sign(dy(i)) > 0
            
            y_new(i) = y(i) + delta_y(i,1);
        else
            y_new(i) = y(i);
        end
    end

    theta_new(i) = mod(atan2(y_new(i), x_new(i)), 2*pi)  ;   
    
    %%% Stop at ID %%%
    
    if sqrt(x_new(i)^2 + y_new(i)^2) >= ri
       
        x_new(i) = ri*cos(theta_new(i)); 
        y_new(i) = ri*sin(theta_new(i)); 
        
    end 
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%   CALCULATE ARC LENGTH  %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

L = 0; 

for i = 2:length(x_new)
    
   L = L + (sqrt((x_new(i) - x_new(i-1))^2 + (y_new(i) - y_new(i-1))^2) + sqrt((x(i) - x(i-1))^2 + (y(i) - y(i-1))^2))/2; %Takes the average of the old profile and the new profile to calculate the arc length
    
   % Need to remove the arc length value if that position has reached
   % the combustion chamber wall
   
   if sqrt(x_new(i)^2 + y_new(i)^2) >= (ri*0.995)
      
       L = L - (sqrt((x_new(i) - x_new(i-1))^2 + (y_new(i) - y_new(i-1))^2) + sqrt((x(i) - x(i-1))^2 + (y(i) - y(i-1))^2))/2;
       
   end
   
end

A = polyarea(x_new, y_new); %Open area (m^2)






