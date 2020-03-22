function [E_Res, E_Tang] = element_assembly(element_iterator)

global n_en %number of displacement element nodes
global n_ee % number of equations of one element
global node_coords % coordinates of the  nodes 
global IEN % element nodes matrix 
global n_quad % number of quad points for pressure elements
% global dim
global hist_old %history variables between time points
global hist_new %history variables between time points
global chi       
global C_m       
global sigma_iso 
global sigma_ani 
global fiber1_dir
global dt
global t_n1

N = zeros(n_ee,1); %vector
% N_u_symgrad      = zeros(2,2,n_ee_u);
% N_u_div          = zeros(n_ee_u,1);

E_Res=  zeros(n_ee,1);
E_Tang= zeros(n_ee,n_ee);

% [E_soln_new, E_soln_old]=get_element_hom_solution(element_iterator);
[E_soln_new, E_soln_old]=get_element_solution(element_iterator);
% [E_BC_new, E_BC_old]=get_element_BC_solution(element_iterator);

E_I_stim= get_nodal_E_I_stim(element_iterator);%I_stim at the nodes of elem
   
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
    sigma_tens= sigma_iso*eye(2) + sigma_ani*(fiber1_dir*fiber1_dir');
    %quad_coords= get_quad_point_coords(quad_iterator, n_quad, ...
    %    n_en,element_coords);
    
    %get I_stim at this quadrature point
    %Force= get_force_vector(mu, quad_coords_u(1), quad_coords_u(2));
    I_stim=0; 
    %get_stimulus();
    dPhi_I_stim=0;
    
    %Call Material routine, to get:
    %electrical potential rhs  : I_ion,
    %its tangent : dPhi_I_ion
    %recovery variable : r, internal variable  (not necessary)
    [I_ion,dPhi_I_ion,r_new] =  ...
    material_routine(Phi_new, r_old, dt);
    % there is a difference in signs between goktepe (material routine)
    % and the kirshnamoorthi  implementations
    I_ion      = -I_ion; 
    dPhi_I_ion = -dPhi_I_ion;
    
    I_m=I_stim-chi*I_ion;
    dPhi_I_m= dPhi_I_stim - chi*dPhi_I_ion; %does chi depend on Phi?

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
                ( chi* C_m* 1/dt * N(i) * N(j) ...
                - dPhi_I_m * N(i) * N(j) ) ...
                + N_grad(i,:) * sigma_tens * N_grad(j,:)' ...
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
            ( chi* C_m* 1/dt* N(i)* (Phi_new-Phi_old) ...
              - N(i)* I_m ...
              + N_grad(i,:)*sigma_tens*Phi_grad_new ) ...
            * JxW ...
            + E_Res(i) ;
    end
    
end
end