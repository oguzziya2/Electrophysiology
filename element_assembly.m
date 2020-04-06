function [E_Res, E_Tang] = element_assembly(element_iterator)

global n_en %number of displacement element nodes
global n_ee % number of equations of one element
global node_coords % coordinates of the  nodes 
global IEN % element nodes matrix 
global n_quad % number of quad points for pressure elements
% global dim
% global hist_old %history variables between time points
% global hist_new %history variables between time points
global chi       
global C_m       
global sigma_iso 
global sigma_ani 
global fiber1_dir
global dt
% global t_n1

N = zeros(n_ee,1); %vector
% N_u_symgrad      = zeros(2,2,n_ee_u);
% N_u_div          = zeros(n_ee_u,1);

E_Res=  zeros(n_ee,1);
E_Tang= zeros(n_ee,n_ee);

% [E_soln_new, E_soln_old]=get_element_hom_solution(element_iterator);
[E_soln_new, E_soln_old]=get_element_solution(element_iterator);
% [E_BC_new, E_BC_old]=get_element_BC_solution(element_iterator);

%get stimulus  applied at the nodes - traction boundary condition
[E_I_stim, E_dPhi_I_stim]= get_nodal_E_I_stim(element_iterator);%I_stim at the nodes of elem
    
%get the ionic currents calculated for each node, in global assembly
[E_I_ion,E_dPhi_I_ion] = get_E_I_ion(element_iterator);

%calculate I_m from I_ion and I_stim
E_I_m=E_I_stim - chi*E_I_ion;
E_dPhi_I_m=E_dPhi_I_stim - chi*E_dPhi_I_ion; %does chi depend on Phi?


nodes_of_element= IEN (:,element_iterator);
element_coords = node_coords(nodes_of_element,:);

%loop on all quadrature points
%WARNING! : some quadrature points are gauss points some are nodal points
%because of lumped mass matrix of ionic currents
%assume there are same number of quadrature points in both methods

for quad_iterator=1:n_quad
    
    %N is a vector for all shape functions, only this quad point
    N=get_shape_fnc_vals(quad_iterator,n_quad, n_en,"gauss");
    N_lum=get_shape_fnc_vals(quad_iterator,n_quad, n_en,"lumped");
    
    %get shape fnc symmetric gradient and divergence
    N_grad= get_shape_fnc_grad(quad_iterator,n_quad,...
        n_en, element_coords, "gauss");
    N_grad_lum= get_shape_fnc_grad(quad_iterator,n_quad,...
        n_en, element_coords,"lumped");
%    N_div= get_shape_fnc_div(quad_iterator,n_quad, n_en, element_coords);
    
    %calculate Phi and its gradient, at quad point
    %dont use this for surface integration    
    Phi_new = N'*E_soln_new;
    Phi_old = N'*E_soln_old; 
    Phi_grad_new=N_grad'*E_soln_new;
    
    %get  jacobian and integration weights, 
    %not to be used for surface integration
    JxW=get_JxW(quad_iterator,n_quad,n_en,element_coords, "gauss");
    JxW_lum=get_JxW(quad_iterator,n_quad,n_en,element_coords, "lumped");

    %calculate I_m at this quadrature point from element nodal values
    I_m_lum=N_lum'*E_I_m; %results in a scalar
    dPhi_I_m_lum=N_lum.*E_dPhi_I_m; %results in a vector
    
    %get conductivity tensor at this quad point
    sigma_tens= sigma_iso*eye(2) + sigma_ani*(fiber1_dir*fiber1_dir');
    %quad_coords= get_quad_point_coords(quad_iterator, n_quad, ...
    %    n_en,element_coords);
    
    %-----------------------
    %assemble  element tangent
    for node_i=1:n_en
        i=node_i;
        
        for node_j=1:n_en
            j =node_j;
 
            E_Tang(i,j) = ...
                ( chi* C_m* 1/dt * N(i) * N(j) ...
                + N_grad(i,:) * sigma_tens * N_grad(j,:)') ...
                * JxW  ...
                - N_lum(i) * dPhi_I_m_lum(j) * JxW_lum ...%lumped ionic cur
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
              + N_grad(i,:)*sigma_tens*Phi_grad_new ) ...
            * JxW ...
            - N_lum(i)* I_m_lum * JxW_lum ...%lumped ionic cur
            + E_Res(i) ;
    end
    
end
end