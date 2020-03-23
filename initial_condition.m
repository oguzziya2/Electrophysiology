function  initial_condition(face, initial_value)
%
global G_soln_n
global node_coords
global ID
global n_eq   %total number of equations of the sysyem

G_soln_n =zeros(n_eq,1)-80 ; %initiate domain from -80mV

switch face 
    case 'none'
        IC_node_nums=[];
    case 'left' % for left face x_coord is zero
        [IC_node_nums, ~]=find(node_coords(:,1)==0);%x_coord is zero
    otherwise
        error('choice of face is not implemented yet')
end
   
    IC_IDs= ID(IC_node_nums);
    G_soln_n(IC_IDs)=initial_value;

end

