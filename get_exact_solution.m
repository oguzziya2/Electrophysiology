function [V] = get_exact_solution(x ,y, t)
%this solution vector is from sec 4.1.1 of 
% " isogeometric analysis : stable elements for the 2d stokes equation"
% by Buffa et al. IJNMF 2011, 65:1407-1422

%its force vector is  defined in get_force_vector function. 

V= 1/4*(1-cos(2*pi*x))*(1-cos(2*pi*y))*t-80;

end

