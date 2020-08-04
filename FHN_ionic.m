function [f_Phi, dp_fp, r_new] = ...
    FHN_ionic(Phi_new, r_old, dt)
%
 global t_n1
fh_1=65;
fh_2=35;
fh_3=200;
alpha = -0.5; 
a= 0.0;
b= -0.6;
c=50;


%non-dimensionalize
Phi_new_nd=(Phi_new+fh_2)/fh_1;
dt_nd= dt/fh_3;

%from the paper- there is  typo
% r_new  = (r_old + (Phi_new_nd+a)*dt_nd)/(1+b*dt_nd) ;
% f_Phi  = c*(Phi_new_nd*(Phi_new_nd-alpha)*(1-Phi_new_nd) - r_new);
% dp_fp  = c*(-3*Phi_new_nd.^2 +2*(1-alpha)*Phi_new_nd-alpha) ...
%                  -(c*dt_nd)/(1+b*dt_nd);
   
%deal ii implementation 
% r_new  = (r_old-dt_nd/c*(Phi_new_nd-a))/(1+b*dt_nd/c);
% f_Phi = c*(r_new +Phi_new_nd - 1/3*Phi_new_nd.^3 + alpha);
% dp_fp = c*(1-Phi_new_nd.^2)-(c*dt_nd)/(c+dt_nd*b);

%from the paper goktepe2009- my hand derivation:
r_new = (dt_nd*(Phi_new_nd+a)+r_old)/(b*dt_nd +1);
f_Phi = c*(-Phi_new_nd^3 + Phi_new_nd^2*(alpha+1) - alpha*Phi_new_nd - r_new);
dp_fp = c*(-3*Phi_new_nd^2 + 2*(alpha+1)*Phi_new_nd -alpha) ...
                                                    -(c*dt_nd)/(b*dt_nd+1);


% 
% %RE- DIMENSIONALIZE
f_Phi=fh_1*f_Phi/fh_3 ;
dp_fp=dp_fp/fh_3;


% %for manuf soln testing with constant fphi
% f_Phi=10;
% dp_fp=0;
% r_new= t_n1;

end

