function probe_node_id= preprocess (input_val, shape, renumbering_type)

global n_ed 
global n_el 
global n_en
% global n_en_u
% global n_en_p
global n_ee
global n_el_x
% global n_el_y
% global n_ee_u
% global n_ee_p
global n_np
global LM
global IEN
global ID
global BC
global node_coords
% global dim

n_ed=1;  %number of element dofs
n_en=4;%number of  element nodes (4, 8 or 9)
probe_coord=[ 2.5, 0.0]; % probe coord should be a node coordinate size 1x2

n_ee=n_en;


if (shape=="square")
    n_el_x=input_val;
    create_square_mesh();
elseif (shape=="rectangle")
    n_el_x=input_val;
    create_rectangular_mesh(2.5,0.1);% rectangular mesh of 2.5x0.1 cm
else
    error('choice of shape is not there yet')
end

if ( exist('renumbering_type','var') ) 
    %renumber global eqns  according to  displacement and pressure comeponents
    index_map=dof_renumbering(renumbering_type);
end


% CREATE LOCATION MATRIX
LM=zeros(n_ed, n_en, n_el);
for i=1:n_ed
    for a=1:n_en
        for e=1:n_el
            LM(i,a,e)= ID(i,IEN(a,e));
        end
    end
end

probe_node_id = ID(find_probe_node(probe_coord));

%-----------------------------------------------
%DOF TO  EQN AND EQN TO DOF ARE NOT USED ANYMORE
%eqn numbers corresponding to  node numbers are hard coded into
% global assembly and element assembly routines, 
%-----------------------------------------------



%CREATE EQN_TO_DOF AND DOF_TO_EQN for 9 node q2q1 element, 3 dofs (u,u,p)
% Q2Q1 9 node velocity- 4 node pressure element.     
%      dofs:(u1,u2,p)
%
%
%  4(7,8,22)    7(13,14)     3(5,6,21)
% o------------o------------o
% |            |            |
% |8(15,16)    |9(17,18)    |6(11,12)
% o------------o------------o
% |            |            |
% |1(1,2,19)   |5(9,10)     |2(3,4,20)
% o------------o------------o

% %use: 
% %[equation  number ]= dof_to_eqn ( dof, node number  )
% dof_to_eqn= [ 1  3  5  7  9  11 13 15 17 ;
%               2  4  6  8  10 12 14 16 18 ;
%               19 20 21 22 0  0  0  0  0   ];
% %use:          
% % [dof, node number, component index ] = eqn_to_dof (equation number , :)          
% eqn_to_dof= [ 1 1 1;
%               2 1 2;
%               1 2 3;
%               2 2 4;
%               1 3 5;
%               2 3 6;
%               1 4 7;
%               2 4 8;
%               1 5 9;
%               2 5 10;
%               1 6 11;
%               2 6 12;
%               1 7 13;
%               2 7 14;
%               1 8 15;
%               2 8 16;
%               1 9 17;
%               2 9 18;
%               3 1 1;
%               3 2 2;
%               3 3 3;
%               3 4 4; ];

end
