function [I_stim] = get_I_stim(x, y, t)

% %return I_stim at the point (x,y) at time t
% if((x<0.1)&&(y<0.1)&&(t<30))
%     I_stim = 1000;
% else
%     I_stim =0 ;
% end

% %below is for FitzHugh Nagumo model manufactured soln
I_stim= (43750*cos(2*pi*x))/9 - (175*t)/12 + (43750*cos(2*pi*y))/9 ...
    - (175*exp((3*t)/1000)*(250*cos(2*pi*x) + 250*cos(2*pi*y) ...
    - 250*cos(2*pi*x)*cos(2*pi*y) - 115))/9 + (175*t*cos(2*pi*x))/12 ... 
    + (175*t*cos(2*pi*y))/12 - (43750*cos(2*pi*x)*cos(2*pi*y))/9 ... 
    + 140*(cos(2*pi*x)/4 - 1/4)*(cos(2*pi*y) - 1) ...
    + 2275*((t*(cos(2*pi*x)/4 - 1/4)*(cos(2*pi*y) - 1))/65 - 9/13) ...
    *((t*(cos(2*pi*x)/4 - 1/4)*(cos(2*pi*y) - 1))/65 - 5/26) ... 
    *((t*(cos(2*pi*x)/4 - 1/4)*(cos(2*pi*y) - 1))/65 - 22/13) ... 
    - (175*t*cos(2*pi*x)*cos(2*pi*y))/12 + (2*t*pi^2*cos(2*pi*y) ... 
    *(cos(2*pi*x)/4 - 1/4))/5 + (t*pi^2*cos(2*pi*x)*(cos(2*pi*y) - 1))/10 - 20125/9;




%below is for constant Iion =10 for manufactured soln
% I_stim= 140*(cos(2*pi*x)/4 - 1/4)*(cos(2*pi*y) - 1) ...
%     + (2*t*pi^2*cos(2*pi*y)*(cos(2*pi*x)/4 - 1/4))/5 ...
%     + (t*pi^2*cos(2*pi*x)*(cos(2*pi*y) - 1))/10 - 1400 ;

end

