function [E_BC_new, E_BC_old]=get_element_BC_solution(element_iterator)

global n_ee
global IEN
global LM
global BC

E_BC_new=zeros(n_ee,1);
E_BC_old=zeros(n_ee,1);

for i=1:size(E_BC_new,1) %local node iterator
    
    j=1;  %dof number for electrical potential
    %!!!!!!!!!!!!!!!!!!!!IMPORTANT TO CHANGE IT AGAIN 
    % IN ELECTROMECHANICS IMPLEMENTATION
    
    global_eqn_index= LM(j,i,element_iterator);
    
    if (global_eqn_index~=0)
        E_BC_new(i,1) =  0 ;
        E_BC_old(i,1) =  0 ;

    elseif (global_eqn_index==0)
        global_node_num=IEN(i, element_iterator);
        E_BC_new(i,1)=BC(global_node_num,j);
        E_BC_old(i,1)=BC(global_node_num,j);

    else
        error('coulnt find solution value for this elemnt node')
    end
    
end

end
