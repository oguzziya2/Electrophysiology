%This code solves electrophysiology problem using quad elements in 2d
% Oguz Ziya Tikenogullari - 2019

clc
clear all
close all


% global n_ed %number of element dofs
global n_el %number of elements
% global n_en % max  of n_en 's
%global n_en_u %number of displacement element nodes
%global n_en_p  %number of pressure element nodes
% global n_ee % number of equations of one element
%global n_ee_u % number of disp eqns for one element
%global n_ee_p % number of pressure eqns for one element
global n_eq   %total number of equations of the sysyem
%global n_quad_p % number of quad points for pressure elements
%global n_quad_u  % number of quad points for disp elements
global n_quad %number of quad points 
% global n_np  %number of nodal points
global n_h %number of history variables per quad point.
% global LM    %Loccation matrix
global G_soln_n %global solution matrix
global G_soln_n1 %global solution matrix
% global nodal_I_stim % global stimulus values at the nodes 
% global IEN % element nodes matrix 
global ID % destination matrix
% global BC % dirichlet boundary  conditions
% global node_coords % coordinates of the  nodes 
global dim  % number of spatial dimensions
global hist_old %history variables between time points
global hist_new %history variables between time points
global chi       
global C_m       
global sigma_iso 
global sigma_ani 
global fiber1_dir
global dt %time step size
global t_n1
global tol

dim           =2;     	     % number of spatial dimensions

%problem parameters 
% probe_node   =             %node  number to check activation
chi           =140;   	     %cm^-1
C_m           =1;      	     %microFarat/cm^2
sigma_iso     =0.1;%0.0176; 	     %mS/cm
sigma_ani     =0;%0.1158; 	     %mS/cm  %TRANSVERSELY ANISOTROPIC 
fiber1_dir    =zeros(dim,1); 
fiber1_dir(1) =1;            %1 fiber family  in 2 dimensions 
%I_stim       =50000; 	     %microA/cc applied for 2 ms 
t_final       =100;          %time to finalize simulation
dt            =1;     	     %time step size-fixed
tol           =1e-8;  	     %tolerance for norm  check of global residual
newton_maxi   =10;    	     %maximum number of newton iterations 
n_quad        =4;            %number of quadrature points
n_h           =1;            %number of internal variables
flag          =false;        %first activation flag

% square ([0,1]x[0,1]) or rectangle ([0, 2.5]x[0,.1])
% input value is the element number in x-dir
probe_node_id= preprocess (20,'square'); 
%probe node is where we check the activation time  at

n_eq = max(ID,[],'all'); %number of global equations 

%initialize variables
t_n      =0;             %t_0
t_n1     =0;  

% set initial conditions
initial_condition('none',-80); 
%none : no initial condition:everywhere is at -80 mV

hist_old = zeros(size(G_soln_n)); %r_0
hist_new = zeros(size(G_soln_n)); %for gp storage zeros(n_el, n_quad, n_h)

%set nodal tractions (I_stim)
set_nodal_I_stim('left',0); %stimulus from left boundary

%initialize solution vector at t_n+1
G_soln_n1=G_soln_n ;

output_results(t_n1);%output initial solution 

%time marching
while (t_n1<t_final-tol)
    
    %new time step 
    t_n =t_n1;
    t_n1=t_n+dt; 
    %fprintf("time step: %d \n", t_n1);
    
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
        
        %fprintf("Residual of norm: %e \n", Norm_Res);

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
    
    % output solution variables at selected time points
    output_results(t_n1);
    
    %check  the activation of probe node 
%     [activated, tmp ]=is_active(probe_node_id);
%     if (activated && (flag==0))
%         activation_time= tmp;
%         flag=1;
%     end
    
    % update G_soln .
    G_soln_n = G_soln_n1 ; 
    
    % update  history variable
    hist_old = hist_new ;
    
end
disp('FIN')

% %calculate errors 
% [L2_error_velocity, L2_error_pressure]= postprocess();
