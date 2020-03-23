function create_rectangular_mesh () 
%this one creates square mesh for Q1 elements 

 
global n_el 
global n_el_x
global n_el_y
global n_np
global n_en
global n_ed
global BC
global IEN
global ID
global node_coords


%CREATE THE MESH, BC, IEN, ID
% n_el_x= 2;
n_el_y= 1;
 
n_el=n_el_x*n_el_y;

n_rows= n_el_y+1;
n_cols= n_el_x+1;

row_coords= linspace(0,0.1,n_rows); %dimensions of rectangular mesh is
col_coords= linspace(0,2.5,n_cols); % 2.5x0.1 cm

% node_coords = zeros( 2, n_rows*n_cols, 'double');
% ID = zeros( 1, n_rows*n_cols, 'double');

for row=1:n_rows
       
    %create nodal coordinates
    add_n_coords= [ col_coords; repelem(row_coords(row),n_cols)];
    node_coords = [ node_coords add_n_coords];
        
    %create ID matrix from ones
    add_ID=zeros(n_ed,n_cols);
    if (row_coords(row) ~= 0) %If bottom row, zeros(dirichlet)
        add_ID(1,:) = 1;
    else 
        add_ID(1,:) = 1;%If bottom row, zeros(dirichlet)
    end
    ID= [ID add_ID] ;
    
    %create IEN 
    add_node_numbers=ones(1, n_cols);    
    node_numbers(row,:) = add_node_numbers;
end 
node_numbers=node_numbers';

% assign equation numbers to ID
[i,j]=find(ID) ; 
for k=1:length(j)
    ID(i(k), j(k)) = k ;
end
% ID hasn't finished yet


%assign node numbers  
[i, j]=find(node_numbers) ;
for k=1:length(j)
    node_numbers(i(k), j(k)) = k ;
end
node_numbers=node_numbers';

[n_rows,n_cols]=size(node_numbers);
k=0;
for i=1:(n_rows-1)
    for j=1:(n_cols-1)
        k=k+1; %  element number counter 
        col_index=[j j+1  j+1  j   ];
        row_index=[i i    i+1  i+1 ];
        for l=1:n_en
            IEN(l,k)=node_numbers(row_index(l), col_index(l));
        end
    end
end

[~,n_np]=size(node_coords); %number of nodal points

BC= zeros (n_np, n_ed);
boundary_idx= find(ID==0);

% fix Pressure to  zero at one node
BC(boundary_idx) = 0; %indices have swithced

node_coords=node_coords'; %we need the transpose for the rest of the code. 
%fix this 
end














