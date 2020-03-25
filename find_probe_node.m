function [probe_node] = find_probe_node(probe_coord)

global node_coords
global dim

if ( (size(probe_coord,1) ~= 1) ||...
        (size(probe_coord,2) ~= size (node_coords,2)))
    error('Probe coord is not right size')
end

tmp = (node_coords==probe_coord);
tmp2 = ( (tmp(:,1)&tmp(:,2)) == 1); 
probe_node =find(tmp2);
end

