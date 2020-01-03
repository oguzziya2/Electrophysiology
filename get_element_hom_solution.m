function [E_soln_new, E_soln_old] = get_element_hom_solution(element_iterator)

global n_ee
%global IEN
global LM
%global BC
global G_soln_n
global G_soln_n1

E_soln_new=zeros(n_ee,1);
E_soln_old=zeros(n_ee,1);

for i=1:size(E_soln_new,1) %local node iterator
    
    j=1;  %dof number for electrical potential
    %!!!!!!!!!!!!!!!!!!!!IMPORTANT TO CHANGE IT AGAIN 
    % IN ELECTROMECHANICS IMPLEMENTATION
    
    global_eqn_index= LM(j,i,element_iterator);
    
    if (global_eqn_index~=0)
        E_soln_new(i,1) = G_soln_n1(global_eqn_index);
        E_soln_old(i,1) = G_soln_n(global_eqn_index);

    elseif (global_eqn_index==0)
        %global_node_num=IEN(i, element_iterator);
        E_soln_new(i,1)=0;
        E_soln_old(i,1)=0;

    else
        error('coulnt find solution value for this elemnt node')
    end
    
end

end
