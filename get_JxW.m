function [JxW] = get_JxW(quad_iterator, n_quad, element_node_number, ...
    node_coordinates, quad_rule)

if ~exist('quad_rule','var')
    % default quad rule is gauss
    error("quad rule does not exist")
end

dx_dxi= get_dx_dxi (quad_iterator, n_quad, ...
    element_node_number, node_coordinates, quad_rule);

J=det(dx_dxi);

if (quad_rule == "lumped")
    
    switch n_quad
        case 4
            w1=1; w2=1;
            weights=[w1*w1 w1*w2 w2*w1 w2*w2];
            W=weights(quad_iterator);
            
        case 9
            error ("check this weights of nodal quad for 9 node") 
            
        otherwise
            error ('this quad rule is not implemented')
    end
   
elseif (quad_rule=="gauss")
    
    switch n_quad
        case 4
            w1=1; w2=1;
            weights=[w1*w1 w1*w2 w2*w1 w2*w2];
            W=weights(quad_iterator);
            
        case 9
            w1=5/9; w2=8/9; w3=5/9;
            weights=[ w1*w1 w1*w2 w1*w3 ...
                w2*w1 w2*w2 w2*w3 ...
                w3*w1 w3*w2 w3*w3    ];
            W=weights(quad_iterator);
            
        otherwise
            error ('this quad rule is not implemented')
    end
end 

JxW= J * W ; 

end

