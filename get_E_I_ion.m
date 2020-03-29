function [E_I_ion,E_dPhi_I_ion] = get_E_I_ion(element_iterator)

% get the ionic current and its derivative for the new time point
global I_ion
global dPhi_I_ion
global n_ee
global LM

E_I_ion=zeros(n_ee,1); 
E_dPhi_I_ion=zeros(n_ee,1); 

%the time dependency of applied stimulus 

for i=1:size(E_I_ion,1) %local node iterator
    
    j=1;  %dof number for electrical potential
    %!!!!!!!!!!!!!!!!!!!!IMPORTANT TO CHANGE IT AGAIN 
    % IN ELECTROMECHANICS IMPLEMENTATION
    
    global_eqn_index= LM(j,i,element_iterator);
    if (global_eqn_index~=0)
        E_I_ion(i,1) = I_ion(global_eqn_index);
        E_dPhi_I_ion(i,1) = dPhi_I_ion(global_eqn_index);

    elseif (global_eqn_index==0) %this is not used in EP b/c no Drichlet
        E_I_ion(i,1) =0 ;
        E_dPhi_I_ion(i,1) =0 ;
        error('there should not be history varable at the dirichlet node')
    else
        error('coulnt find  value for this elemnt node')
    end
end

end

