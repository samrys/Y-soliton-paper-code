function KP_solver_single( t, Lx, Nx, Nt,...
                                Ly, Ny,...
                                u0,...
                                data_dir )
% Solves: KP eq. (u_t + uu_x + u_xxx)_x + u_yy = 0
% on [-xmax,xmax] & [-ymax,ymax] by FFT in space with
% Window technique described in Kao, Kodama 2010 in y-direction
% Pseudospectral in space
% integrating factor v = exp[+i(k^3*epsilon^2-l^2/k)t]*u_hat 
% and RK4 in time
%
%
% Inputs:
%
% t        :  1D array of output times
% Lx, Nx   :  spatial grid parameters (x)
% Ly, Ny   :  spatial grid parameters (y)
% u0     :  structure containing IC and asymptotic approximation 
%             of solution and its derivs over time
% data_dir : output directory where data at each output timestep
%             is written to a file #####.mat
%
% Outputs:  NONE except for data written to files

%% Set up temporal vector, increment counter, data directory
global tout inc dir
    tout = single(t);
    inc  = 0;
    dir  = data_dir;
    dt   = single(10^-Nt);

%% Define spatial and wavenumber grids
domain = struct;
   domain.x = single((2*Lx/Nx)*(-Nx/2:Nx/2-1)');
   domain.dx = single((2*Lx/Nx));
   domain.Lx = single(Lx);
   domain.y = single((2*Ly/Ny)*(-Ny/2:Ny/2-1)');
   domain.dy = single((2*Ly/Ny));
   domain.Ly = single(Ly);
   [domain.X,domain.Y] = meshgrid(domain.x,domain.y);
   domain.k = single((pi/Lx)*[0:Nx/2-1 0 -Nx/2+1:-1]');
   domain.l = single((pi/Ly)*[0:Ny/2-1 0 -Ny/2+1:-1]');
   [domain.KX,domain.KY] = meshgrid(domain.k,domain.l);
   domain.invKx = 1./(1i*domain.KX);
    
%% Write static matrices here, so not constantly redefining  
    %% Windowing function (from Kao 2010), rescaled, and derivs
    n = 27; an = single((1.111)^n*log(10));
    W.o = single(exp( -an * abs(domain.Y/Ly).^n ));
    W.p = single(-an*n/Ly * W.o .* abs(domain.Y/Ly).^(n-1).*sign(domain.Y/Ly));
    W.pp = single(an*n/Ly^2 * W.o .* abs(domain.Y/Ly).^(n-2) .* ...
               ( (-(n-1) + an*n*abs(domain.Y/Ly).^n).*sign(domain.Y/Ly).^2 ));
    iphi = zeros(size(domain.KX),'single');
    iphi(domain.KX~=0) = 1i*(domain.KY(domain.KX~=0).^2./(domain.KX(domain.KX~=0)) - domain.KX(domain.KX~=0).^3);
    
    %% Integrating factor matrics also precomputed
    domain.E0 = exp(iphi*tout(1));
    domain.Em0 = exp(-iphi*tout(1));
    domain.E1 = exp(0.5*dt*iphi);
    domain.Em1 = exp(-0.5*dt*iphi);

    
%% Save variables, prepping to call the solver
    % Output what we are about to do
    disp(['Solving KP eqtn.']);
    disp([' Time interval:  [', num2str(t(1)),...
           ',',num2str(t(end)),']']);
    % Construct initial condition on spatial domain
    u_init = u0.u0;
    u      = u_init; tnow = min(t);
    % Construct initial v, the windowed version of u
    v_init = single(u_init - (1 - W.o) .* u0.u0);
    v      = v_init;
    % Construct initial Vhat, the fft'd, integrating factored v
    Vhat_init = fft2(v_init);
    Vhat      = Vhat_init;
    % Hardcode in what Vhat should be for KX = 0
%     domain.VhK0        = Vhat(domain.KX==0);
    
    % Save initial condition (inc = 0)
	save(strcat(data_dir,num2str(inc,'%05d')),...
                'u','u_init','v','v_init','Vhat',...
                'Vhat_init','inc','tnow');
	inc = inc + 1;
    disp(['Time = ',num2str(tout(inc)),', inc = ',int2str(inc),...
          '/',int2str(length(tout))]);
	start = tic;

%% Call the solver
    myRK4_single( Vhat_init, u0, dt, W,...
                    domain,start)
 
%% Finish and clean up
    finish = toc(start);
    disp('Calculation Finished');
    days_left = datenum([0 0 0 0 0 toc(start)]);
    time_left=datevec(days_left-floor(days_left));
    disp(['Computation time = ',...
        int2str(floor(days_left)),'d ',...
        int2str(time_left(4)),'h ',...
        int2str(time_left(5)),'m ',...
        num2str(time_left(6)),'s']);
save([data_dir,'parameters.mat'],'finish','-append'); 