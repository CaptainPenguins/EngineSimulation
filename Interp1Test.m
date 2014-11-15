x = 0:pi/4:10*pi;
v = sin(x);
xq = 0:pi/16:2*pi;
x = mod(x, 2*pi);

vq1 = interp1(x,v,xq, 'spline');


figure
plot(x,v,'o',xq,vq1,':.');
xlim([0 4*pi]);
title('(Default) Linear Interpolation');