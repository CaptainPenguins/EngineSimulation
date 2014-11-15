function [x,y] = CoreShapeGenerator(ri, R, N, IR, OR)
% clear all; close all; 
% 
% ri = 2.5235/25.4; 
% R = 0.5/25.4; %Spoke radius
% N = 5; %Number of spokes
% IR = 1.15/25.4; %Radius of inner circle 
% OR = 1.5/25.4; %Radius of spoke arm

s = 0.0002; %Spacing

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%    MAKE ERROR MESSAGES  %%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (OR+R)>ri
    error('Spoke outer diameter exceeds combustion chamber inner diameter. Reduce spoke outer diameter to less than %0.1f. Aborting', (ri - R)*0.95); 
end

if (IR >= OR)
    error('Inner circle diameter exceeds spoke arm diameter. Reduce inner circle diameter or increase spoke arm length. Aborting'); 
end

if (R <= (0.204/2)/25.4*1.1)
    error('Spoke radius is too small to allow fasteners, increase spoke radius. Aborting'); 
end

if(asin(3*R/2/(R/2+IR)))*2*N > 2*pi
    error('Spoke arms are overlapping with each other. Decrease number of spoke arms. Aborting.'); 
end

if N==1
    error('Number of spoke arms shall not be equal to 1. Aborting.'); 
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%    CREATE SPOKE TIPS  %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Initialize Variables %%%

phiST = zeros(floor(pi/(s/R)), N);
xST = zeros(size(phiST)); yST = zeros(size(phiST));
theta = zeros(N, 1);

for n = 1:N
    
    theta(n) = (n-1)*2*pi/N; %Locates the center of the spokes
    
    phiST(:,n) = (3*pi/2+(n-1)*2*pi/N)+s/(2*R):s/R:(5*pi/2+(n-1)*2*pi/N)-s/(2*R);
    
    xST(:,n) = R*cos(phiST(:,n)) + OR*cos((n-1)*2*pi/N);
    yST(:,n) = R*sin(phiST(:,n)) + OR*sin((n-1)*2*pi/N);
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%   CREATE RADII   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

r = R/2;
phir_L = zeros(ceil((pi/2 - asin((R+r)/(r+IR)))/(s/r)),N);
phir_R = zeros(ceil((pi/2 - asin((R+r)/(r+IR)))/(s/r)),N);

xr_L = zeros(size(phir_L)); yr_L = zeros(size(phir_L));
xr_R = zeros(size(phir_L)); yr_R = zeros(size(phir_L));
xr = zeros(2*size(phir_L,1), N); yr = zeros(2*size(phir_L,1), N);

for n = 1:N
    
    phir_L(:,n) = (pi + theta(n) + asin((R+r)/(r+IR))) : s/r : 3*pi/2 + theta(n);
    phir_R(:,n) = (pi/2 + theta(n)) : s/r : (pi - asin((R+r)/(r+IR)) + theta(n));
    
    xr_L(:,n) = r*cos(phir_L(:,n)) + (r+IR)*cos(asin((R+r)/(r+IR)) + theta(n));
    yr_L(:,n) = r*sin(phir_L(:,n)) + (r+IR)*sin(asin((R+r)/(r+IR)) + theta(n));
    
    xr_R(:,n) = r*cos(phir_R(:,n)) + (r+IR)*cos(-asin((R+r)/(r+IR)) + theta(n));
    yr_R(:,n) = r*sin(phir_R(:,n)) + (r+IR)*sin(-asin((R+r)/(r+IR)) + theta(n));
    
    xr(:,n) = [xr_L(:,n); xr_R(:,n)];
    yr(:,n) = [yr_L(:,n); yr_R(:,n)];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%   CREATE LINES   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

xl_L = zeros(floor((xST(end,1) - xr_L(end,1))/s),N); yl_L = zeros(size(xl_L));
xl_R = zeros(size(xl_L)); yl_R = zeros(size(yl_L));

xl = zeros(length(xl_L)*2, N);
yl = zeros(size(xl));

for n = 1:N
    
    if abs(sin(theta(n))) <= 1e-10
        xl_L(:,n) = xr_L(end,n) + s*cos(theta(n))/2: s*cos(theta(n)) : xST(end,n) - s*cos(theta(n))/2;
        yl_L(:,n) = yr_L(end,n): (yST(end,n) - yr_L(end,n))/(length(xl_L(:,n))-1) : yST(end,n) ;
    else if abs(cos(theta(n))) <= 1e-10
            yl_L(:,n) = yr_L(end,n) + s*sin(theta(n))/2: s*sin(theta(n)) : yST(end,n) - s*sin(theta(n))/2;
            xl_L(:,n) = xr_L(end,n): (xST(end,n) - xr_L(end,n))/(length(yl_L(:,n))-1) : xST(end,n) ;
            
        else
            xl_L(:,n) = xr_L(end,n) + s*cos(theta(n))/2: s*cos(theta(n)) : xST(end,n) - s*cos(theta(n))/2;
            yl_L(:,n) = yr_L(end,n) + s*sin(theta(n))/2: s*sin(theta(n)) : yST(end,n) - s*sin(theta(n))/2;
        end
    end
    
    xl_R(:,n) = xl_L(:,n) + 2*R*sin(theta(n));
    yl_R(:,n) = yl_L(:,n) - 2*R*cos(theta(n));
    
    xl(:,n) = [xl_L(:,n);xl_R(:,n)];
    yl(:,n) = [yl_L(:,n);yl_R(:,n)];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%   CREATE CIRCLE   %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

phiC = zeros(ceil((theta(2) - 2*asin((R+r)/(r+IR)))/(s/IR)), N);

xC = zeros(size(phiC)); yC = zeros(size(phiC));

for n = 1:N
    
    if n == N
        phiC(:,n) = asin((R+r)/(r+IR)) + theta(n): s/IR : 2*pi - asin((R+r)/(r+IR));
    else
        phiC(:,n) = asin((R+r)/(r+IR)) + theta(n): s/IR : theta(n+1) - asin((R+r)/(r+IR));
    end
    
    xC(:,n) = IR*cos(phiC(:,n));
    yC(:,n) = IR*sin(phiC(:,n));
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%   COMBINE AND SORT ARRAYS   %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

xUS = [xST(:);xr(:);xl(:);xC(:)];
yUS = [yST(:);yr(:);yl(:);yC(:)];

phi = mod(atan2(yUS,xUS), 2*pi);
r = sqrt(xUS.^2+yUS.^2);

A = [r, phi];
B = sortrows(A, 2);

u = unique(B(:,2)); 
n = hist(B(:,2), u); 
p = u(n>1);


for i = 1:length(p)
   
    indices = find(B(:,2) == p(i)); 
    B(indices(1),:) = []; 
    
end

x = zeros(length(B(:,2)), 1); 
y = zeros(size(x)); 

for zeta = 1:length(B(:,2))    
    x(zeta) = B(zeta,1) * cos(B(zeta, 2));
    y(zeta) = B(zeta,1) * sin(B(zeta, 2));
end

%%% Plot Wall for Reference %%%

phi_wall = 0:0.01:2*pi; 
x_wall = ri*cos(phi_wall); 
y_wall = ri*sin(phi_wall); 


plot(x,y, 'r', x_wall, y_wall, 'b'); axis equal; legend('Core Shape', 'Chamber Wall'); 


