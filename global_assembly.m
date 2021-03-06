function [G_Res, G_Tang] = global_assembly()

global mu
global n_ed %number of element dofs
global n_el %number of elements
global n_en % max  of n_en 's
global n_ee % number of equations of one element
global n_eq
global n_np  %number of nodal points
global LM    %Loccation matrix
global G_soln_n %global solution matrix
global G_soln_n1 %global solution matrix
global IEN % element nodes matrix 
global ID % destination matrix
global BC % dirichlet boundary  conditions
global node_coords % coordinates of the  nodes 
global dim  % number of spatial dimensions

% n_eq=max(ID(:));  % number of global equations

G_Res            = zeros(n_eq,1);
G_Tang           = zeros(n_eq,n_eq);

for element_iterator=1:n_el
    
    %get element right hand side and element tangent matrices
    [E_Res, E_Tang] = element_assembly(element_iterator);
    %[E_Res, E_Tang] = test_assembly(n_ee,0);
    
    %assemble element tangent and residual to global ones
    % i,j are eqn indices in element level, of the node of concern
    % I,J are  "    "      " global  level, of the node of concern
    
    %iterate over the nodes of the element
    for node_iterator_i=1:n_en
        
        %one equation per node for electrophysiology,
        % so , i = node_iterator
        i=node_iterator_i ;
        I=LM(n_ed,node_iterator_i,element_iterator);
        
        if (I==0)% if this dof corresponds to dirichlet BC, skip
            continue
        end
        
        for node_iterator_j=1:n_en
            
            %one equation per node for electrophysiology,
            % so , i = node_iterator
            j=node_iterator_j ;
            J=LM(n_ed,node_iterator_j,element_iterator);
            
            if (J==0) %if this dof corresponds to dirichlet BC, skip
                continue
            end
            
            G_Tang(I,J) = G_Tang(I,J) ...
                + E_Tang(i,j);
        end
        G_Res(I) = G_Res(I)+ E_Res(i);
    end
end

end

