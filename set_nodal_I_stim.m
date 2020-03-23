function set_nodal_I_stim(face, stim_value)

%set the stimulus position in this function,
%set time dependency from the get_nodal_E_I_stim function.

global nodal_I_stim %global stimulus values at the nodes
global n_eq   %total number of equations of the system
global node_coords % coordinates of the  nodes 
global ID % destination matrix

% n_nodes= size(ID,2); %number of nodes of the fe mesh
% nodal_I_stim= zeros(n_eq,1);

switch face 
    case 'left' % for left face x_coord <= 0.1
        [stim_node_nums, ~]=find(node_coords(:,1)<=0.1);%x_coord <= 0.1 
    otherwise
        error('choice of face is not implemented yet')
end


stim_IDs= ID(stim_node_nums);
nodal_I_stim=zeros(n_eq,1);
nodal_I_stim(stim_IDs)= stim_value;


end

