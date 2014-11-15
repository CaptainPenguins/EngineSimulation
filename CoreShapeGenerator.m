function [x,y] = CoreShapeGenerator(ri, R_s, N_s, r_ic, r_s, C)


%% Author: J. Osborne, October 2013

% This code will generate the (x,y) coordinates - in metres -  for either a wagon-wheel
% or a circular fuel core given a set of input parameters. The input parameters are:

%ri = Engine inside radius (m)
%R_s = Half-width of spoke (m)
%N_s = number of spokes 
%r_ic = Radius of inner core (m) 
%r_s = Radial distance to center of spoke (m)
%C = Geometry generation selector (0 = wagon wheel, 1 = circle)

if C == 0
    
    
    s = 0.0002; %Spacing
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%    MAKE ERROR MESSAGES  %%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if (r_s+R_s)>ri
        error('Spoke outer diameter exceeds combustion chamber inner diameter. Reduce spoke outer diameter to less than %0.1f. Aborting', (ri - R_s)*0.95);
    end
    
    if (r_ic >= r_s)
        error('Inner circle diameter exceeds spoke arm diameter. Reduce inner circle diameter or increase spoke arm length. Aborting');
    end
    
    if (R_s <= (0.204/2)/25.4*1.1)
        error('Spoke radius is too small, increase spoke radius. Aborting');
    end
    
    if(asin(3*R_s/2/(R_s/2+r_ic)))*2*N_s > 2*pi
        error('Spoke arms are overlapping with each other. Decrease number of spoke arms. Aborting.');
    end
    
    if N_s==1
        error('Number of spoke arms shall not be equal to 1. Aborting.');
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%    CREATE SPOKE TIPS  %%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%% Initialize Variables %%%
    
    phiST = zeros(floor(pi/(s/R_s)), N_s);
    xST = zeros(size(phiST)); yST = zeros(size(phiST));
    theta = zeros(N_s, 1);
    
    for n = 1:N_s
        
        theta(n) = (n-1)*2*pi/N_s; %Locates the center of the spokes
        
        phiST(:,n) = (3*pi/2+(n-1)*2*pi/N_s)+s/(2*R_s):s/R_s:(5*pi/2+(n-1)*2*pi/N_s)-s/(2*R_s);
        
        xST(:,n) = R_s*cos(phiST(:,n)) + r_s*cos((n-1)*2*pi/N_s);
        yST(:,n) = R_s*sin(phiST(:,n)) + r_s*sin((n-1)*2*pi/N_s);
        
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%   CREATE RADII   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    r = R_s/2;
    phir_L = zeros(ceil((pi/2 - asin((R_s+r)/(r+r_ic)))/(s/r)),N_s);
    phir_R = zeros(ceil((pi/2 - asin((R_s+r)/(r+r_ic)))/(s/r)),N_s);
    
    xr_L = zeros(size(phir_L)); yr_L = zeros(size(phir_L));
    xr_R = zeros(size(phir_L)); yr_R = zeros(size(phir_L));
    xr = zeros(2*size(phir_L,1), N_s); yr = zeros(2*size(phir_L,1), N_s);
    
    for n = 1:N_s
        
        phir_L(:,n) = (pi + theta(n) + asin((R_s+r)/(r+r_ic))) : s/r : 3*pi/2 + theta(n);
        phir_R(:,n) = (pi/2 + theta(n)) : s/r : (pi - asin((R_s+r)/(r+r_ic)) + theta(n));
        
        xr_L(:,n) = r*cos(phir_L(:,n)) + (r+r_ic)*cos(asin((R_s+r)/(r+r_ic)) + theta(n));
        yr_L(:,n) = r*sin(phir_L(:,n)) + (r+r_ic)*sin(asin((R_s+r)/(r+r_ic)) + theta(n));
        
        xr_R(:,n) = r*cos(phir_R(:,n)) + (r+r_ic)*cos(-asin((R_s+r)/(r+r_ic)) + theta(n));
        yr_R(:,n) = r*sin(phir_R(:,n)) + (r+r_ic)*sin(-asin((R_s+r)/(r+r_ic)) + theta(n));
        
        xr(:,n) = [xr_L(:,n); xr_R(:,n)];
        yr(:,n) = [yr_L(:,n); yr_R(:,n)];
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%   CREATE LINES   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    xl_L = zeros(floor((xST(end,1) - xr_L(end,1))/s),N_s); yl_L = zeros(size(xl_L));
    xl_R = zeros(size(xl_L)); yl_R = zeros(size(yl_L));
    
    xl = zeros(length(xl_L)*2, N_s);
    yl = zeros(size(xl));
    
    for n = 1:N_s
        
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
        
        xl_R(:,n) = xl_L(:,n) + 2*R_s*sin(theta(n));
        yl_R(:,n) = yl_L(:,n) - 2*R_s*cos(theta(n));
        
        xl(:,n) = [xl_L(:,n);xl_R(:,n)];
        yl(:,n) = [yl_L(:,n);yl_R(:,n)];
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%   CREATE CIRCLE   %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    phiC = zeros(ceil((theta(2) - 2*asin((R_s+r)/(r+r_ic)))/(s/r_ic)), N_s);
    
    xC = zeros(size(phiC)); yC = zeros(size(phiC));
    
    for n = 1:N_s
        
        if n == N_s
            phiC(:,n) = asin((R_s+r)/(r+r_ic)) + theta(n): s/r_ic : 2*pi - asin((R_s+r)/(r+r_ic));
        else
            phiC(:,n) = asin((R_s+r)/(r+r_ic)) + theta(n): s/r_ic : theta(n+1) - asin((R_s+r)/(r+r_ic));
        end
        
        xC(:,n) = r_ic*cos(phiC(:,n));
        yC(:,n) = r_ic*sin(phiC(:,n));
        
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
    
else if C == 1
        
        theta = 0:0.01:2*pi;
        x = r_s * sin(theta);
        y = r_s * cos(theta);
    end
end




