function [E_Res, E_Tang] = element_assembly(element_iterator)

global n_en %number of displacement element nodes
global n_ee % number of equations of one element
global node_coords % coordinates of the  nodes 
global IEN % element nodes matrix 
global n_quad % number of quad points for pressure elements
global dim
global hist_old %history variables between time points
global hist_new %history variables between time points
global d_iso   %conductivity coefficient
global dt

N = zeros(n_ee,1); %vector
% N_u_symgrad      = zeros(2,2,n_ee_u);
% N_u_div          = zeros(n_ee_u,1);

E_Res=  zeros(n_ee,1);
E_Tang= zeros(n_ee,n_ee);

% [E_soln_new, E_soln_old]=get_element_hom_solution(element_iterator);
[E_soln_new, E_soln_old]=get_element_solution(element_iterator);
% [E_BC_new, E_BC_old]=get_element_BC_solution(element_iterator);
nodes_of_element= IEN (:,element_iterator);
element_coords = node_coords(nodes_of_element,:);

%calculate residual of pressure
%loop on all quadrature points
for quad_iterator=1:n_quad
    
    %get the history variable (r_old)
    r_old = hist_old(element_iterator, quad_iterator);
    
    %N is a vector for all shape functions, only this quad point
    N=get_shape_fnc_vals(quad_iterator,n_quad, n_en);
   
    %get shape fnc symmetric gradient and divergence
    N_grad= get_shape_fnc_grad(quad_iterator,n_quad,...
        n_en, element_coords);
%     N_div= get_shape_fnc_div(quad_iterator,n_quad, n_en, element_coords);
    
    %calculate Phi and its gradient, at quad point
    %dont use this at the right hand side    
    Phi_new = N'*E_soln_new;
    Phi_old = N'*E_soln_old; 
    Phi_grad_new=N_grad'*E_soln_new;
    
    %get  jacobian and integration weights
    JxW=get_JxW(quad_iterator,n_quad,n_en,element_coords);
    
    %get conductivity tensor at this quad point
    D_tens= d_iso*eye(2);
    %quad_coords= get_quad_point_coords(quad_iterator, n_quad, ...
    %    n_en,element_coords);
    
    %get body force vector at this quadrature point
    %Force= get_force_vector(mu, quad_coords_u(1), quad_coords_u(2));
    
    %Call Material routine, to get:
    %electrical potential rhs  : f_Phi,
    %its tangent : dp_fp
    %recovery variable : internal variable  (not necessary)
    [f_p,dp_fp,r_new] =  ...
    material_routine(Phi_new, r_old, dt);

    %update history variables
    hist_new(element_iterator, quad_iterator,1) = r_new;
    %-----------------------
 
    
    %-----------------------
    %assemble  element tangent
    for node_i=1:n_en
        i=node_i;
        
        for node_j=1:n_en
            j =node_j;
 
            E_Tang(i,j) = ...
                ( 1/dt * N(i) * N(j) ...
                + N_grad(i,:) * D_tens * N_grad(j,:)' ...
                - dp_fp * N(i) * N(j) ) ...
                *JxW  ...
                +E_Tang(i,j) ;
        end
    end
    
     
    %-----------------------
  
    
    %-----------------------
    %assemble  element right hand side
    for node_i=1:n_en
        i=node_i;
        E_Res(i)= ...
            ( N(i)* ((Phi_new-Phi_old)/ dt - f_p) ...
              + N_grad(i,:)*D_tens*Phi_grad_new ) ...
            * JxW ...
            + E_Res(i) ;
    end
    
end
end