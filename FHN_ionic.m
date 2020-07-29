function [f_Phi, dp_fp, r_new] = ...
    FHN_ionic(Phi_new, r_old, dt)
%
% global t_n1
fh_1=35;
fh_2=30;
fh_3=20;
alpha = -0.4; 
a= 0.7;
b= 0.8;
c= 3;


%non-dimensionalize
Phi_new_nd=(Phi_new+fh_2)/fh_1;
dt_nd= dt/fh_3;

% r_new  = (r_old + (Phi_new_nd+a)*dt_nd)/(1+b*dt_nd) ;
r_new  = (r_old-dt_nd/c*(Phi_new_nd-a))/(1+b*dt_nd/c);
% f_Phi  = c*(Phi_new_nd*(Phi_new_nd-alpha)*(1-Phi_new_nd) - r_new);
f_Phi = c*(r_new +Phi_new_nd - 1/3*Phi_new_nd.^3 + alpha);
% dp_fp  = c*(-3*Phi_new_nd.^2 +2*(1-alpha)*Phi_new_nd-alpha) ...
%                  -(c*dt_nd)/(1+b*dt_nd);
dp_fp = c*(1-Phi_new_nd.^2)-(c*dt_nd)/(c+dt_nd*b);
             
% 
% %RE- DIMENSIONALIZE
f_Phi=fh_1*f_Phi/fh_3 ;
dp_fp=dp_fp/fh_3;

end

