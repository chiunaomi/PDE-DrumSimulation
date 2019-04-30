clear all
clf

A = 10000;
alpha = .01;
rho = 0.05; 
rstar = 0.5; 
thetastar = 0;

x = 0:0.1:20;
J = zeros(5,201);
for i = 0:4
    J(i+1,:) = besselj(i,x);
end
bessellambda = find_zeros(J);

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

t = 1; 
T = exp(-alpha*t/2).*(cos((sqrt(alpha^2+4.*bessellambda).*t)./2)+sin((sqrt(alpha^2+4.*bessellambda).*t)./2));
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

