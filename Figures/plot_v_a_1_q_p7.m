% Plot results from Mach reflection soliton simulation
data_dir = 'C:\Users\samry\Documents\MATLAB\KP2\data\miles_resonant\v_shape_a_1_q_p7\';
load([data_dir,'\parameters.mat'],'Lx','Ly');

ran = [0 2.9]; %Range for caxis, placed here for convenience

% Set up spatial grid
fac = 1; % For Fourier interpolation
Nx = fac*Lx*4; Ny = fac*Ly*4;
x0 = 0; % Initial position of soli segment
y0 = 0;
dx = 2*Lx/Nx;
x = dx*[-Nx/2:Nx/2-1];
dy = 2*Ly/Ny;
y = dy*[-Ny/2:Ny/2-1];
[foo,xmin] = min(abs(x-(x0-20)));
[foo,xmax] = min(abs(x-(x0+420)));
[foo,ymin] = min(abs(y-(y0-245)));
[foo,ymax] = min(abs(y-(y0+245)));
xp = unique(round(linspace(xmin,xmax,250)));
yp = unique(round(linspace(ymin,ymax,250)));

% Set up temporal grid
dt = 10; 
tout = [0,100,250]/dt;
Np = length(tout);

% Setup figure
ml = 0.13; % Margin left
mr = 0.03; % Margin right
mt = 0.04; % Margin top
mb = 0.24;  % Margin bottom
pb = 0.01; % Interaxes padding bottom
pr = 0.00; % Interaxes padding right
spanx = (1-ml-mr-(Np-1)*pr)/Np;
spany = (1-mt-mb);
fig_width = 16; % in cm
fig_height = 6;
fontsize = 14;
cmap = load('CoolWarmFloat257.csv');

fh=figure(1);
clf();
fh.Renderer = 'Painters';
set(gcf,'Resize','off')
set(fh,'paperposition',[0,0,fig_width,fig_height],...
       'papersize',[fig_width,fig_height],'paperunits',...
       'centimeters','units','centimeters');

for ii=Np
    inc = tout(ii);
    load([data_dir,'/',sprintf('%05d',inc),'.mat'],'u');

    %Interpolate if desired
    if fac > 1
        u = fftInterpolate(u,fac*size(u));
    end
    
    axes('Position',[ml,mb,spanx,spany]);
    contourf(x(xp)-x0,y(yp),u(yp,xp),25,'edgecolor','none');
    colormap(cmap);
    caxis(ran);
    set(gca,'FontSize',fontsize,'TickLabelInterpreter','latex',...
            'ytick',round([-200,0,200]),...
            'xtick',[-200,0,200,400]);
    xlabel('$x$','interpreter','latex');
    if ii == Np
        ylabel('$y$','interpreter','latex');
    else
        set(gca,'yticklabel',{'','',''});
    end
    axis([x(xmin)-x0,round(x(xmax)-x0),round(y(ymin)),round(y(ymax))]);
end

set(fh,'paperposition',[0,0,fig_width,fig_height],...
       'papersize',[fig_width,fig_height],'paperunits',...
       'centimeters','units','centimeters');

 doc_name = ['v_shape_a_1_q_p7_',num2str(tout(1)*dt),...
     '_',num2str(tout(2)*dt),'_',num2str(tout(3)*dt),'.pdf'];  
 print(fh,'-dpdf',doc_name);


