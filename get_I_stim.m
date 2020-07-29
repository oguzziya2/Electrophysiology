function [I_stim] = get_I_stim(x, y, t)
%return I_stim at the point (x,y) at time t
% if((x<0.1)&&(y<0.1)&&(t<30))
%     I_stim = 1000;
% else 
%     I_stim =0 ;
% end

%below is for constant Iion =10 for manufactured soln
I_stim= 140*(cos(2*pi*x)/4 - 1/4)*(cos(2*pi*y) - 1) ...
    + (2*t*pi^2*cos(2*pi*y)*(cos(2*pi*x)/4 - 1/4))/5 ...
    + (t*pi^2*cos(2*pi*x)*(cos(2*pi*y) - 1))/10 - 1400 ;

end

