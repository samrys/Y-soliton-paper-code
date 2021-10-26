%% Driver file for KP2 numerical solver:  kink initial conditions
%% Domain, ICs, etc go here
save_on  = 1;  % Set to nonzero if you want to run the solver, set
               % to 0 if you want to plot
rmdir_on = 1;  % Set to nonzero if you want to delete and remake the chosen directory
               % Useful for debugging
mail_on  = 1;  % set to nonzero to send email when completed
               % set to 0 to not run email                
plot_on  = 0;  % Set to 1 if you want to plot just before and just
               % after (possibly) calling the solver          
check_IC = 0;  % Set to nonzero to plot the ICs and BCs without running the solver



    %% Numerical Parameters
    tmax   = 200;      % Solver will run from t=0 to t = tmax
    numout = tmax/10+1; % numout times will be saved (including ICs)
    Lx     = 1024;     % Solver will run on x \in [-Lx,Lx]
    Ly     = 512;     % Solver will run on y \in [-Ly,Ly]
    Nx     = Lx*4;    % Number of Fourier modes in x-direction
    Ny     = Ly*4;    % Number of Fourier modes in y-direction

    t      = linspace(0,tmax,numout);
    Nt      = 3;
    dt      = 10^(-Nt); %Timestep size
    s       = 1;   %Size of transition in tanh for q(y)
 
    %% Bent Soliton parameters
    sa = 1;
    qu = -1.4;
    ql = -qu;
    x0 = 0;
    
    xplot  = (2*Lx/Nx)*[-Nx/2:Nx/2-1];
    yplot  = (2*Ly/Ny)*[-Ny/2:Ny/2-1];
    [XPLOT,YPLOT] = meshgrid(xplot,yplot);
      
    qmat = (qu-ql)/2*tanh(YPLOT./s) + (ql+qu)/2;          
    qymat = (qu-ql)/2*1/s*sech(1/s*(YPLOT)).^2;
    intqy = s*(qu-ql)/2*(log(cosh(YPLOT/s)))+YPLOT.*(ql+qu)/2;
    soli = kink(XPLOT,YPLOT,sa,qmat,qymat,intqy,x0);
     
     %Preallocate Calculations
 
        ic_type = ['_soli_kink_','_sa_',num2str(sa),...
                '_qu_',num2str(qu),'_ql_',num2str(ql),...
                '_s_',num2str(s)];

    %% Generate directory, save parameters
	q = strsplit(pwd,filesep);
    %% determine data storage location
    if strcmp(computer,'MACI64')
        maindir = '/Volumes/Data Storage/Numerics/KP';
        if ~exist(maindir,'dir')
            q = strsplit(pwd,filesep);
            maindir = strjoin(q(1:end-1),filesep);
        end
    elseif strcmp(computer,'PCWIN64')
        maindir = 'C:\Users\samry\Documents\MATLAB\KP2';
    else
        q = strsplit(pwd,filesep);
        maindir = strjoin(q(1:end-1),filesep);
        path(path,[strjoin(q(1:end-1),filesep),filesep,'KP_Scratch',filesep,'Solver'])
    end
    
    slant = filesep;


    %% Create directory run will be saved to
    data_dir = [maindir,slant,'Numerics',slant,'KP',slant,...
                '_kink_',...
                '_tmax_',   num2str(round(tmax)),...
                '_Lx_',     num2str(Lx),...
                '_Nx_',     num2str(Nx),...
                '_Ly_',     num2str(Ly),...
                '_Ny_',     num2str(Ny),...
                '_x0_',     num2str(x0),...
                '_init_condns_',ic_type,...
                slant];
    % Create the data directory if necessary
    if ~exist(data_dir,'dir')
        mkdir(data_dir);
    else
        disp(['Warning, directory ',data_dir]);
        if rmdir_on == 1
            disp('already exists, rewriting entire folder');
            rmdir(data_dir,'s')
            mkdir(data_dir);
        else
            disp('already exists, possibly overwriting data');
        end
    end

    savefile = sprintf('%sparameters.mat',data_dir);

    %% If chosen, run the solver using the parameters and conditions above
    if save_on
        if plot_on
            % Load initial data
            
              tplot  = linspace(0,tmax,floor(tmax*10));
              theta = soli.th(0);
              u_init = soli.ua(theta);
              uax = soli.uax(theta).*u_init;
              uay = soli.uay(0).*uax;

            % Plot initial conditions and boundary conditions
            fontsize = 12;
            figure(1); clf;
            subplot(2,2,1)
                contourf(XPLOT,YPLOT,u_init,100,'edgecolor','none'); xlabel('x'); ylabel('y'); 
                title('Initial Conditions');
            subplot(2,2,2)
                contourf(XPLOT,YPLOT,u_init,100,'edgecolor','none'); xlabel('x'); ylabel('y'); 
                title('Asymptotic u');
            subplot(2,2,3)
                contourf(XPLOT,YPLOT,uax,100,'edgecolor','none'); xlabel('x'); ylabel('y'); 
                title('Asymptotic u, x-deriv');
            subplot(2,2,4)
                contourf(XPLOT,YPLOT,uay,100,'edgecolor','none'); xlabel('x'); ylabel('y'); 
                title('Asymptotic u, y-deriv');
            set(gca,'fontsize',fontsize,'fontname','times');
            pause(0.25);
            if check_IC
                drawnow; 
                 return;
            end
        end

        % Save parameters
            save(savefile,'t','Nx','Lx',...
                              'Ny','Ly','Nt');
        % Run timestepper
            KP_solver_single( t, Lx, Nx, Nt,...
                                   Ly, Ny,...
                                   soli,...
                                   data_dir );     
    else
        load(savefile);
    end


if plot_on
    plot_data_fun_2D(data_dir);
    figure(4);
    print('sim','-dpng');

end

if mail_on
      send_mail_message('samuel.ryskamp','Matlab',['Simulation ',data_dir,'done'])
end     
 %% Function that helps construct the initial conditions and the asymptotic solution
    % The values need to be multiplied correctly to be accurate
    % Code ensures minimal time spent evaluating be reusing values
    
function [ soli ] = kink(x,y, sa, qmat, qymat,intqy,x0)
    soli.th = @(t) sa/sqrt(12).*((x-x0) + intqy - (sa^2/3+qmat.^2) .* t);                          
    soli.ua = @(z) sa^2*(sech(z)).^2;
    soli.uax = @(z) -sa/sqrt(3).*tanh(z);
    soli.uay = @(t) sqrt(12)*(y.*qymat+qmat-2.*qmat.*qymat.*t)*(sa/sqrt(12));
    th0 = soli.th(0);
    soli.u0 = soli.ua(th0);
 end     