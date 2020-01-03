function output_results(time_step)
% load('matlab.mat')

global ID 
global BC
global G_soln_n1
global node_coords
global n_el

n_node_x= sqrt(n_el)+1; %number of nodes in  x direction
n_node_y= sqrt(n_el)+1; %number of nodes in  y direction

%2 dimensions, so Z component is zero
X=reshape(node_coords(:,1), [n_node_x, n_node_y ] ) ;
Y=reshape(node_coords(:,2), [n_node_x, n_node_y ] ) ;

Z=zeros(size(X));
r=zeros(size(node_coords,1),1);

for i=1:length (ID)
    if (ID(i) ~= 0)
        r(i)= G_soln_n1(ID(i));
    else 
        r(i)=BC(i);
    end
    
end
% r=[ 1 2 3 4 5 6 7 8 9] ;
% s = mesh(X,Y,Z,'FaceAlpha','0.5')

filename = sprintf('Solution/SOL_%d.vtk',time_step);
vtkwrite(filename,'structured_grid',X,Y,Z,'scalars','title',r);

end
