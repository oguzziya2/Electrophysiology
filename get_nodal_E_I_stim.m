function [E_I_stim,E_dPhi_I_stim] = get_nodal_E_I_stim(element_iterator)

% get the I_stim values for the element pointed by the element iterator

global t_n1         % current time value
global nodal_I_stim %global stimulus values at the nodes
global n_ee
global LM

E_I_stim=zeros(n_ee,1);
E_dPhi_I_stim=zeros(n_ee,1); %this stays as zero because it's a dead load

%the time dependency of applied stimulus 

%apply first 5 ms
if (t_n1 <= 5 ) 
for i=1:size(E_I_stim,1) %local node iterator
    
    j=1;  %dof number for electrical potential
    %!!!!!!!!!!!!!!!!!!!!IMPORTANT TO CHANGE IT AGAIN 
    % IN ELECTROMECHANICS IMPLEMENTATION
    
    global_eqn_index= LM(j,i,element_iterator);
    if (global_eqn_index~=0)
        E_I_stim(i,1) = nodal_I_stim(global_eqn_index);
    elseif (global_eqn_index==0) 
        E_I_stim(i,1) =0 ;
        error('this is not used in EP b/c no Drichlet')
    else
        error('coulnt find stimulus value for this elemnt node')
    end
end
end

end

