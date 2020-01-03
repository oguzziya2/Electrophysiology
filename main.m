%This code solves electrophysiology problem using quad elements in 2d
% Oguz Ziya Tikenogullari - 2019

clc
clear all
close all


% global mu
global n_ed %number of element dofs
global n_el %number of elements
global n_en % max  of n_en 's
%global n_en_u %number of displacement element nodes
%global n_en_p  %number of pressure element nodes
global n_ee % number of equations of one element
%global n_ee_u % number of disp eqns for one element
%global n_ee_p % number of pressure eqns for one element
global n_eq
%global n_quad_p % number of quad points for pressure elements
%global n_quad_u  % number of quad points for disp elements
global n_quad %number of quad points 
global n_np  %number of nodal points
global n_h %number of history variables per quad point.
global LM    %Loccation matrix
global G_soln_n %global solution matrix
global G_soln_n1 %global solution matrix
global IEN % element nodes matrix 
global ID % destination matrix
global BC % dirichlet boundary  conditions
global node_coords % coordinates of the  nodes 
global dim  % number of spatial dimensions
global hist_old %history variables between time points
global hist_new %history variables between time points
global d_iso    %isotropic component of conductivity tensor
global dt %time step size
global tol

%problem parameters 
%mu = 0.04; 
dim      =2;     % number of spatial dimensions
t_final  =500;     %time to finalize simulation
dt       =1;   %time step size-fixed
tol      =1e-8;  %tolerance for norm  check of global residual
newton_maxi=10;  %maximum number of newton iterations 
n_quad   = 4; 
n_h      = 1;
d_iso    = 0.001; 

% square [0,1]x[0,1]
%0: 1 element, 
%1: 4 element, 
%2: 16 element,
%3: 64 element,
%4: 256 element,
%5: 1024 element,
%6: 1 element all dofs free
%7: 1 element test case, clamped 4 sides 
preprocess (2); 
n_eq = max(ID,[],'all'); %number of global equations 


%initialize variables
t_n      =0;             %t_0
t_n1     =0;             
G_soln_n =zeros(n_eq,1)-80 ; %Phi_0
G_soln_n(1:(sqrt(n_el)+1),1) = 0; %bottom edge is at zero initial cond.
G_soln_n1=G_soln_n ;
hist_old = zeros(n_el, n_quad, n_h); %r_0
hist_new = zeros(n_el, n_quad, n_h);

output_results(t_n1);%output initial solution 

%time marching
while (t_n1<t_final-tol)
    
    %new time step 
    t_n =t_n1;
    t_n1=t_n+dt; 
    
    G_soln_n  = G_soln_n1;
    %G_soln_n1= ?? aim is to  find this in this time increment
        
    %global newton iterations, because of nonlinear nature of the problem
    for newton_iter = 0:newton_maxi
        
        %assemble the  global  tangent and residual
        [G_Res, G_Tang] = global_assembly();
        %[G_Res, G_Tang] = test_assembly (n_eq, newton_iter);
        
        % Check for convergence with newly updated residual
        % break newton loop if error is small
        Norm_Res= norm(G_Res,2);
        if(Norm_Res < tol)
            break
        end
        
        % %if solution hasn't converged yet:
        % %solve the system to find newton increment of solution vector
        % G_Soln is global variable.
        % Soln_inc = global_solve(G_Res, G_Tang);
        Soln_inc= -G_Tang\G_Res; 
        
        %update G_soln_n1
        G_soln_n1 = G_soln_n1 + Soln_inc ; 
        
    end
    
    % optional: output solution variables at selected time points
    output_results(t_n1);
    % update G_soln .
    G_soln_n = G_soln_n1 ; 
    
    % update  history variable
    hist_old = hist_new ;
    
end

% %calculate errors 
% [L2_error_velocity, L2_error_pressure]= postprocess();
