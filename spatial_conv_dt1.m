%elements_in_x-dir    activation_time
%
%
%25		         130.4
%50		         134.4
%100		     135.83
%200		     136.45
%400		     136.7478
%800		     136.8462

mesh_size =1./[25 50 100 200 400 800];

%-----------ionic cur at node, GP integration -----------%
node_gp_conduction_vel= 2.5 ./ ...
    ([130.3415 134.2835 135.8100 136.4477 136.7642 136.8958] ./1000);

node_gp_exact_vel = 2.5 ./ (136.8958./1000);
node_gp_error= abs(node_gp_conduction_vel-node_gp_exact_vel)/node_gp_exact_vel;

%-----------classic gauss point-----------%
gp_conduction_vel= 2.5 ./ ([130.4 134.4 135.83 136.45 136.7478 136.8462] ./1000);

gp_exact_vel = 2.5 ./ (136.8462./1000);
gp_error= abs(gp_conduction_vel-gp_exact_vel)/gp_exact_vel;

%-----------ionic cur at node, nodal integr => lumped mass-----------%
lumped_conduction_vel= 2.5 ./ ([144.2235 137.8885 136.6876 136.6672 136.8032 136.8032] ./1000);

lumped_exact_vel = 2.5 ./ (136.8032./1000);
lumped_error= abs(lumped_conduction_vel-lumped_exact_vel)/lumped_exact_vel;

%-----------higher order integr, (GP+NP)/2-----------%
higher_conduction_vel= 2.5 ./ ([136.5393 135.8209 136.1931 136.5416 136.7708 136.9021] ./1000);

higher_exact_vel = 2.5 ./ (136.9021./1000);
higher_error= abs(higher_conduction_vel-higher_exact_vel)/higher_exact_vel;
%===========PLOT===========%
figure(1)
loglog(mesh_size, node_gp_error,'b-o',...
       mesh_size, gp_error, 'r-o',    ...
       mesh_size, lumped_error, 'g-o', ...
       mesh_size, higher_error, 'c-o');
xlabel('mesh size [mm]')
ylabel('relative error of conduction velocity')
legend('currents at node, GP integr.',...
       'currents at Gauss point' ,...
       'currents at node, nodal integ.', ...
       'currents at node, higher order integ.')

figure(2)
semilogx( mesh_size, node_gp_conduction_vel, 'b-o' , ...
          mesh_size, gp_conduction_vel, 'r-o' ,      ...
          mesh_size, lumped_conduction_vel, 'g-o', ...
          mesh_size, higher_conduction_vel, 'c-o');
xlabel('mesh size [mm]')
ylabel('conduction velocity')
legend('currents at node, GP integr.',...
       'currents at Gauss point' ,...
       'currents at node, nodal integ.', ...
       'currents at node, higher order integ.')