function [pRF,alpha0] = pRF_elliptical(x,y,l_rho,radius,l_theta,l_sigma_x,l_sigma_y,l_alpha,r_min)

% functions definition
rho=@(l_rho,radius) radius*normcdf(l_rho,0,1);

theta=@(l_theta) 2*pi.*normcdf(l_theta,0,1)-pi; 

sigma=@(l_sigma,radius, r_min)(radius-r_min).*normcdf(l_sigma,0,1)+r_min;

miu_x0=@(l_rho,radius,l_theta) rho(l_rho,radius).*cos(theta(l_theta));
       
miu_y0=@(l_rho,radius,l_theta) rho(l_rho,radius).*sin(theta(l_theta));

alpha=@(l_alpha) 2*pi/2.*normcdf(l_alpha,0,1)-pi/2; 

%% elliptic pRF with rotation

sigma_x=sigma(l_sigma_x,radius,r_min);
sigma_y=sigma(l_sigma_y,radius,r_min);
alpha0=alpha(l_alpha);
x0=miu_x0(l_rho, radius, l_theta);
y0=miu_y0(l_rho, radius, l_theta);
a = cos(alpha(l_alpha))^2/(2*sigma_x^2) + sin(alpha(l_alpha))^2/(2*sigma_y^2);
b = -sin(2*alpha(l_alpha))/(4*sigma_x^2) + sin(2*alpha(l_alpha))/(4*sigma_y^2);
c = sin(alpha(l_alpha))^2/(2*sigma_x^2) + cos(alpha(l_alpha))^2/(2*sigma_y^2);
pRF= exp( - (a*(x-x0).^2 - 2*b*(x-x0).*(y-y0) + c*(y-y0).^2));
%figure, imagesc(pRF)



end