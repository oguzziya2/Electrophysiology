function [I_stim] = get_I_stim(x, y, t)
%return I_stim at the point (x,y) at time t
I_stim = 140*(cos(2*pi*x)/4 - 1/4)*(cos(2*pi*y) - 1) ...
    + 56*t*pi^2*cos(2*pi*y)*(cos(2*pi*x)/4 - 1/4) ...
    + 14*t*pi^2*cos(2*pi*x)*(cos(2*pi*y) - 1);
end

