alpha = 10000;
rho = 0.05; 
rstar = 0.5; 
thetastar = 0;

[r,theta] = meshgird(0:0.1:1,0:pi/30;2*pi);
t = 1; 
T = exp(-alpha*t/2)*(cos(sqrt(alpha^2+4*lambda)/2*t)+sin(sqrt(alpha^2+4*lambda)/2*t));
