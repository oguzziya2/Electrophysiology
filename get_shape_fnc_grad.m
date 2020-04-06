function [N_grad] = get_shape_fnc_grad(quad_iterator, n_quad, ...
    element_node_number, node_coordinates, quad_rule)

if ~exist('quad_rule','var')
    % default quad rule is gauss
    error("quad rule does not exist")
end

N_derv= get_shape_fnc_derv(quad_iterator,n_quad, ...
    element_node_number, quad_rule);

dx_dxi= get_dx_dxi(quad_iterator, n_quad, element_node_number, ...
    node_coordinates, quad_rule);

% dxi_dx= inv(dx_dxi);
% N_grad=N_derv*dxi_dx;

N_grad=N_derv/dx_dxi;
end