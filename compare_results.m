function  [error]= compare_results()

global n_el %number of elements
global n_en %number of element nodes
% global n_en_p  %number of pressure element nodes
global dim  % number of spatial dimensions
% global n_ee_u % number of disp eqns for one element
global IEN % element nodes matrix 
global node_coords % coordinates of the  nodes 
global n_quad
global t_n1



N = zeros(n_en,1); %vector
error = 0;

for element_iterator=1:n_el
     
    [E_soln, ~]=get_element_solution(element_iterator);
    
    nodes_of_element= IEN (:,element_iterator);
    element_coords = node_coords(nodes_of_element,:);

    for quad_iterator=1:n_quad

        N=get_shape_fnc_vals(quad_iterator,n_quad, n_en);
        JxW=get_JxW(quad_iterator,n_quad,n_en,element_coords);
        
        quad_coords= get_quad_point_coords(quad_iterator,n_quad, ...
            n_en,element_coords);
        
        V_exact = get_exact_solution(quad_coords(1),quad_coords(2), t_n1);
        
        V= E_soln'* N;
        
        V_diff =V-V_exact;
        
        error= (V_diff^2)*JxW + error;
    end
end

disp('L2 norm of Potential :')
error=sqrt(error)
end

