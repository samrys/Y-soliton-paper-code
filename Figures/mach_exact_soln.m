function [aout1,aout2,qout1,qout2] = mach_exact_soln(yvec,t,a0,q0)
%Exact solution for Mach reflection from modulation theory
aout1 = zeros(size(yvec));
aout2 = aout1;
qout1 = aout1;
qout2 = aout1; 

s = 12;

q2 = sqrt(a0);
a2 = q0^2;
ct = 2/3*(sqrt(a0)-q0)*t;
cs = (2*q2-2/3*sqrt(a2))*t-s;
cz = 2*(q2+sqrt(a2))*t-s;

for ii = 1:length(yvec)
    y = yvec(ii);
    if y>ct
        aout1(ii) = a0;
        qout1(ii) = -q0;       
        if y<cs
            aout2(ii) = a2;
            qout2(ii) = q2;
        elseif y>=cs && y<cz
            aout2(ii) = 9/64*(2*(q2+sqrt(a2))-(y+s)/t)^2;
            qout2(ii) = 1/8*(3*((y+s)/t)+2*(q2+sqrt(a2)));
        else
            aout2(ii) = 0;
            qout2(ii) = q2+sqrt(a2);
        end
    elseif y<-ct
        aout1(ii) = a0;
        qout1(ii) = q0;       
        if y>-cs
            aout2(ii) = a2;
            qout2(ii) = -q2;
        elseif y<=-cs && y>-cz
            aout2(ii) = 9/64*(2*(q2+sqrt(a2))+(y-s)/t)^2;
            qout2(ii) = -1/8*(-3*(y-s)/t+2*(q2+sqrt(a2)));
        else
            aout2(ii) = 0;
            qout2(ii) = -(q2+sqrt(a2));
        end
    else
        aout1(ii) = (q0+sqrt(a0))^2;
        qout1(ii) = 0;
        aout2(ii) = NaN;
        qout2(ii) = NaN;
    end
end
    
