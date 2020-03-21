function  initial_condition(face, initial_value)
%
global G_soln_n
global node_coords
global ID


switch face 
    case 'left'
        [IC_node_nums, ~]=find(node_coords(:,1)==0);%x_coord is zero
    otherwise
        error('choice of mesh is not there yet')
end
   
    IC_IDs= ID(IC_node_nums);
    G_soln_n(IC_IDs)=0;

end

