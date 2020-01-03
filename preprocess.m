function preprocess (mesh_choice , renumbering_type)

global n_ed 
global n_el 
global n_en
% global n_en_u
% global n_en_p
global n_ee
% global n_ee_u
% global n_ee_p
global n_np
global LM
global IEN
global ID
global BC
global node_coords
global dim


n_ed=1;  %number of element dofs
n_en=4;%number of  element nodes (4, 8 or 9)

n_ee=n_en;

switch mesh_choice
    
    case 0 % 1 element test, clamped sides
        create_square_mesh(1);
        
    case 1 % 2x2=4 elements, clamped sides
        create_square_mesh(2);
         
    case 2 %4x4=16 elements, clamped sides
        create_square_mesh(4);
        
    case 3 %8x8=64 elements, clamped sides
        create_square_mesh(8);
        
    case 4 %16x16=256 elements, clamped sides
        create_square_mesh(16);
        
    case 5 %32x32=1024 elements, clamped sides
        create_square_mesh(32);
        
    case 6 %case for test, clamped  bottom 

        n_el=1;  %number of elements
        n_np=4; %number of nodal points
        
        node_coords = [ 0.00 0.00 1.00 1.00 ;
                        0.00 1.00 0.00 1.00 ]' ;
        
        ID = [ 1 2 3 4 ];
        
        IEN=[ 1 3 4 2 ]' ;
        
        BC = zeros(n_np, n_ed);
%         BC(1:2, n_ed)= ??;
        %boundary conditions are all zero
        % to modify, enough to fill non zero boundary conditions
        
%     case 7 %case for test, 1 element clamped 4 sides
%         
%         n_el=1;  %number of elements
%         n_np=9; %number of nodal points
%         
%         node_coords = [ 0.00  0.50 1.00 0.00 0.50 1.00 0.00 0.50 1.00;
%                         0.00  0.00 0.00 0.50 0.50 0.50 1.00 1.00 1.00 ]' ;
%         
%         ID = [0  0  0  0  2  0  0  0  0 ;
%               0  0  0  0  3  0  0  0  0 ;
%               0  0  1  0  0  0  4  0  5 ];
%         
%         IEN=[ 1 3 9 7 2 6 8 4 5 ]' ;
%         
%         BC = zeros(n_np, n_ed);
%         BC(1, 3)= 0.052;
%         %boundary conditions are all zero
%         % to modify, enough to fill non zero boundary conditions
%     
    otherwise
        error('choice of mesh is not there yet')
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
