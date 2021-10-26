function [aout,qout] = reg_exact_soln(yvec,t,a0,q0)

%Exact solution for regular reflection from modulation theory
sh = 0;

qstar = q0+sqrt(a0);
vt = 2*qstar;
vb = 2*q0-2/3*sqrt(a0);



aout = zeros(size(yvec));
qout = zeros(size(yvec));

for ii = 1:length(yvec)
    y = yvec(ii);
    s = sign(y)*sh;
    if abs(y)<= 1
        aout(ii) = 2*a0^2;
        qout(ii) = 0;
    elseif abs(y-s)/t > vt
        aout(ii) = 0;
        qout(ii) = sign(y)*qstar;
    elseif abs(y-s)/t < vb
        aout(ii) = a0;
        qout(ii) = sign(y)*q0;
    else 
        aout(ii) = 9/64*(2*qstar-(sign(y)*(y-s))/t)^2;
        qout(ii) = sign(y)/8*(3*((sign(y)*(y-s))/t)+2*qstar);
    end
end