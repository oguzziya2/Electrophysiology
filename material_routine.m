function [f_Phi,dp_fp,r_new] =  ...
    material_routine(Phi_new, r_old, dt)
%
ap_1=100;
ap_2=80;
ap_3=12.9;
c=8;
alpha =0.01; 

Phi_new_nd=(Phi_new+ap_2)/ap_1;
dt_nd= dt/ap_3;

[r_new, dPhi_r]= compute_r_new(Phi_new_nd, r_old, dt_nd);

f_Phi= c*Phi_new_nd*(Phi_new_nd-alpha) * (1-Phi_new_nd)-r_new*Phi_new_nd;

dp_fp= c*(-3*Phi_new_nd.^2 +2*(1+alpha)*Phi_new_nd-alpha) ...
                -r_new-Phi_new_nd*(dPhi_r);

%RE- DIMENSIONALIZE
f_Phi=ap_1*f_Phi/ap_3 ;
dp_fp=dp_fp/ap_3;

end

function [r_new, dPhi_r]= compute_r_new(Phi_new, r_old ,dt)

  
%constants:
gamma=0.002;
m1=0.2;
m2=0.3;
c=8;
b=0.15;
global tol


r_new= r_old; %initialize
% dPhi_r = 0;
Rr= r_new-r_old-((gamma+(m1*r_new)/(m2+Phi_new)) ...
                      *(-r_new-c*Phi_new*(Phi_new-b-1)))*dt;
%find r_new
icount=0;
while(abs(Rr)>tol)
  icount=icount+1;
  if (icount >10)
    error ("newton loop is not converging")
    break
  end
  dr_Rr=1+(gamma+m1/(m2+Phi_new)*(2*r_new+c*Phi_new*(Phi_new-b-1)))*dt;
  r_new=r_new-Rr/dr_Rr; %update r_new
  Rr= r_new-r_old-((gamma+(m1*r_new)/(m2+Phi_new)) ...
  		         *(-r_new-c*Phi_new*(Phi_new-b-1)))*dt;
end
dr_Rr=1+(gamma+m1/(m2+Phi_new)*(2*r_new+c*Phi_new*(Phi_new-b-1)))*dt;
dp_Rr=((gamma+m1*r_new/(m2+Phi_new))*c*(2*Phi_new-b-1) ...
       -m1*r_new/((m2+Phi_new).^2) ...
       *(r_new+c*Phi_new*(Phi_new-b-1)))*dt;
dPhi_r=-dp_Rr/dr_Rr;        

end
