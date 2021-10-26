function myRK4_single( Vhat_init, u0,...
                    dt, W, domain,abst0)
%myRK4_KP2 is a modified fourth-order explicit Runge-Kutta timestepping method
% Includes variables and saving methods specific to a KP2 solver
% INPUTS:
%   u_init:      Initial condition
%   v_init:      Initial condition, windowed
%   Vhat_init:   Initial condition, windowed, FFT'd, and int-factored
%   soli     :   structure containing IC and asymptotic approximation 
%                   of solution and its derivs over time
%   dt:          Time stepping increment (approximate)
%   tout:        Times to be output
%   W, Wp, Wpp:  Windowing function and its derivatives (exact)
%   iphi:        Integrating factor exponent
%   domain:      Structure containing the following
    %    x,  y:     real space domain (vectors)
    %    X,  Y:     real space domain (matrices)
    %   dx, dy:     Real space discretization in x and y
    %   Lx, Ly:     Real, unscaled space maxima
    %    k,  l:     wavenumber domain (vectors)
    %   KX, KY:     wavenumber domain (matrices)
%   abst0:       time for the beginning of the run!
% OUTPUTS:
%   (none except to file)

    % Call global variables; initialize RK4 method
    global tout inc dir
    Vold = Vhat_init; % fft'd, unshifted, decayed IC
    Ep = domain.E0;
    Em = domain.Em0;

    % Subsequent time steps
    for jj = 1:length(tout)-1
        abstnow = toc(abst0);
        disp(['Calculating ',num2str(jj),' out of ',num2str(length(tout)-1),...
              ' Time: ',...
                int2str(floor(abstnow/60/60)),'h ',...
                int2str(floor(abstnow/60)-floor(abstnow/60/60)*60),'m ',...
                num2str(mod(abstnow,60)),'s']);
            tmid = single(linspace(tout(jj),tout(jj+1),ceil((tout(jj+1)-tout(jj))/dt)+1));
            for ii = 2:length(tmid)
                [Epnew,Emnew,Vnew] = RK4(tmid(ii-1), tmid(ii)-tmid(ii-1), Vold, ...
                            u0,W, domain, Ep, Em );
                        Ep = Epnew;
                        Em = Emnew;
                if sum(isnan(Vnew(:)))>0
                    error(['Not a Number encountered at t=',num2str(tmid(ii))]);
                end

                Vold = Vnew;
            end
        disp('');

        %% Save data
        v    = Vnew.*Em;
        tnow = gather(tout(jj+1));
        th = u0.th(tnow);

        u = gather(ifft2(v,'symmetric') + (1-W.o).*u0.ua(th));
        umax = max(u(ceil(end/2),:));
          save(strcat(dir,num2str(inc,'%05d')),'u','umax','tnow','inc');
          inc = inc +1;
    end


% KP2 RK4 function
function [Epnew,Emnew,Vhatnew] = RK4( t, dt, Vhat, u0, W, domain,Ep,Em )
% Solves: KP eq. (u_t + uu_x + epsilon^2 u_xxx)_x +u_yy = 0
% on [-xmax,xmax] & [-ymax,ymax] by FFT in space with integrating factor 
% v = exp[+i(k^3*epsilon^2-l^2/k)t]*u_hat
%     g  = -1/2*1i*dt*KX;

	E1 = domain.E1;     Em1 = domain.Em1;
    Ezero = Ep;         Ezeroi = Em;
    Ehalf = Ep.*E1;     Ehalfi = Em.*Em1;
    Eone  = Ehalf.*E1;  Eonei  = Ehalfi.*Em1;
    Epnew = Eone;       Emnew = Eonei;

% Function evals of asymptotic solution

        [uazero,uaxzero,uayzero] = u(t,u0);
        [uahalf,uaxhalf,uayhalf] = u(t+dt/2,u0);
        [uaone,uaxone,uayone] = u(t+dt,u0);
        
    Va  = G( Ezero, Ezeroi.* Vhat          ,...
             uazero, uaxzero, uayzero     ,...
             W, domain );
    Vb  = G( Ehalf, Ehalfi.*(Vhat+dt/2*Va) ,...
             uahalf, uaxhalf, uayhalf,...
             W, domain );     % 4th-order
    Vc  = G( Ehalf, Ehalfi.*(Vhat+dt/2*Vb) ,...
             uahalf, uaxhalf, uayhalf,...
             W, domain );     % Runge-Kutta
    Vd  = G( Eone, Eonei  .*(Vhat+dt*Vc) ,...
             uaone  , uaxone  , uayone  ,...
             W, domain );
    Vhatnew = Vhat + dt*(Va + 2*(Vb+Vc) + Vd)/6;

    
% Integrating factor evaluation

function GV = G( Et, vhat, uasy, dxuasy, dyuasy, W, domain )
    v = ifft2(vhat);
    RHS = (1-W.o) .*...
          ifft( 1i*domain.KX.*fft( W.o.*uasy.* dxuasy -...
                 ifft(1i*domain.KX.*fft(v.*uasy,[],2),[],2) ,[],2),[],2)+ ...
          ( 2*W.p.*dyuasy + W.pp.*uasy ) ;
    RHShat = fft(RHS,[],2);
	v2hat  = fft(v.*v,[],2);
    GV  = Et.*( fft(domain.invKx.*RHShat -...
                    1i*domain.KX./2.*v2hat ));
    GV(domain.KX==0) = 0;
    disp('');
    

% Evaluates asymptotic solution at boundaries
 
function [ua, uax, uay] = u(t,u0)
              theta = u0.th(t);
              ua = u0.ua(theta);
              uax = u0.uax(theta).*ua;
              uay = u0.uay(t).*uax;
