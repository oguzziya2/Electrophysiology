clc 
clear all
close all

global tol 
tol = 1e-7;
tmax=  500; 
dt= 1;
t0 =0 ;

time= t0:dt:tmax + dt;

p0=-50;
r0=0;
pold=p0;
rold=r0;

diso=0.1;

phi= zeros (size(time));
rec= zeros (size(time));

time_count =1;

while (time(time_count)<tmax) 
    time_count=time_count+1;
    
    pnew= pold;
    [fp,dp_fp, rnew] =  material_routine(pnew, rold, dt);
    Rp=pnew-pold-dt*fp;
    while (abs(Rp)>tol)
        dp_Rp=1-dt*dp_fp;
        pnew=pnew-Rp/dp_Rp;
        [fp,dp_fp, rnew] =  material_routine(pnew, rold, dt);
        Rp=pnew-pold-dt*fp;
    end
    pold=pnew;
    rold=rnew;
    phi(time_count)=pnew;
    rec(time_count)=rnew;

end

plot(time, phi, time, rec)