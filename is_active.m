function [bool_active, activation_time] = is_active(probe_node_id)
global G_soln_n
global G_soln_n1
global t_n1
global dt

threshold= -40; %activation threshold is -40 mV;

bool_active= G_soln_n1(probe_node_id) >= threshold;

if (bool_active)
    tmp= dt*(G_soln_n1(probe_node_id)-threshold) ...
        /(G_soln_n1(probe_node_id)-G_soln_n(probe_node_id));
    
    activation_time=t_n1-tmp;
else
    activation_time=0;
end


end

