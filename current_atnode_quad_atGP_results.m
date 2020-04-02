%elements_in_x-dir    activation_time
%
%
%25		     130.4
%50		     134.4
%100		     135.83
%200		     136.45
%400		     136.7478
%800		     136.8462

mesh_size =1./[25 50 100 200 400 800];

node_gp_conduction_vel= 2.5 ./ ...
    ([130.3415 134.2835 135.8100 136.4477 136.7642 136.8958] ./1000);

node_gp_exact_vel = 2.5 ./ (136.8958./1000);
node_gp_error= abs(node_gp_conduction_vel-node_gp_exact_vel)/node_gp_exact_vel;

gp_conduction_vel= 2.5 ./ ([130.4 134.4 135.83 136.45 136.7478 136.8462] ./1000);

gp_exact_vel = 2.5 ./ (136.8462./1000);
gp_error= abs(gp_conduction_vel-gp_exact_vel)/gp_exact_vel;

figure(1)
loglog(mesh_size, node_gp_error,'b-o', mesh_size, gp_error, 'r-o');
xlabel('mesh size [mm]')
ylabel('relative error of conduction velocity')
legend('currents at node','currents at Gauss point')

figure(2)
semilogx( mesh_size, node_gp_conduction_vel, 'b-o' , ...
    mesh_size, gp_conduction_vel, 'r-o' )
xlabel('mesh size [mm]')
ylabel('conduction velocity')
legend('currents at node','currents at Gauss point')
