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
conduction_vel= 2.5 ./ ([130.4 134.4 135.83 136.45 136.7478 136.8462] ./1000);

exact_vel = 2.5 ./ (136.8462./1000);
error= abs(conduction_vel-exact_vel)/exact_vel;
loglog(mesh_size, error, '-o');
