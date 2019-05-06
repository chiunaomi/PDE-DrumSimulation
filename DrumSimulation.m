clear all
close all

A = 10000;
alpha = .01;
rho = 0.05; 
rstar = 0.5; 
thetastar = 0;

X = 0:0.1:20;
J = zeros(5,201);
for i = 0:4
    J(i+1,:) = besselj(i,X);
end
bessellambda = find_zeros(J);
for i = 1:5
    if i == 1
        squished(i,:) = X./bessellambda(i,1);
    else
        squished(i,:) = X./bessellambda(i,2);
    end
end
[x,y] = meshgrid(-1:0.1:1,-1:0.1:1);
[theta,r] = cart2pol(x,y);
%r = 0:0.1:1;
%theta = 0:0.2*pi:2*pi;
u = (A/2*pi*rho).*exp((-1/(2*rho^2)).*((r.*cos(theta)-rstar*cos(thetastar)).^2+(r.*sin(theta)-rstar*sin(thetastar)).^2));
[x,y,z] = pol2cart(theta,r,u);
%x = r.*cos(theta);
%y = r.*sin(theta);
figure
% polarplot(x,y);
% hold on
% contourf(x,y,u);
% hold off
surf(x,y,z)

bessellambda = find_zeros(J);
u2 = besselj(10, bessellambda*r)*cos(10*theta);
[x,y,z2] = pol2cart(theta,r,u2);
%x = r.*cos(theta);
%y = r.*sin(theta);
figure
% polarplot(x,y);
% hold on
% contourf(x,y,u);
% hold off
surf(x,y,z2)
t = 1; 
T = exp(-alpha*t/2).*(cos((sqrt(alpha^2+4.*bessellambda).*t)./2)+sin((sqrt(alpha^2+4.*bessellambda).*t)./2));
a = [];
b = [];
for n = 1:10
    base = @(r,theta) cos(n*theta)*besselj(bessellambda(i)*r, 10);
    a(:,i) = integrate(integrate(u(theta, r).*base(r, theta).*r, 0, 1), -pi, pi);
end
    
%% Functions
function lambda = find_zeros(bessel_functions)
    [m,n] = size(bessel_functions);
    lambda = zeros(m,5);
    for i = 1:m
         ind = 1;
        for j = 1:length(bessel_functions)-1
            if ~isequal(sign(bessel_functions(i,j)),sign(bessel_functions(i,j+1)))
                lambda(i,ind) = j*0.1-0.1;
                ind = ind+1;
            end
        end
    end
end

