function [N_div] = get_shape_fnc_div(quad_iterator, n_quad,...
                          element_node_number, node_coordinates, quad_rule)
         
%this is only for div of vector valued field variable, velocity/displacement     

global dim
global n_ee_u

if ~exist('quad_rule','var')
    % default quad rule is gauss
    error("quad rule does not exist")
end

%N_grad= zeros(element_node_number, dim); 

N_grad= get_shape_fnc_grad(quad_iterator,n_quad, element_node_number, ...
                      node_coordinates, quad_rule);
                  
N_div=zeros(n_ee_u,1);

for  i=1:element_node_number
    
    if(dim ~=2)
        error('check following implementation')
    end
    
    N_div(i*dim-1)= N_grad(i,1);
    N_div(i*dim)= N_grad(i,2);
    
end

end

