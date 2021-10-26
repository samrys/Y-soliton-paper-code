% Compare regular reflection simulation results with modulation theory

%Load data from regular reflection simulation
data_dir = 'C:\Users\samry\Documents\MATLAB\KP2\data\miles_resonant\v_shape_a_1_q_1p4';
load([data_dir,'\parameters.mat'],'Lx','Ly');

ran = [0 2]; %Range for caxis, placed here for convenience
a0 = 1;
q0 = 1.4;

% Set up spatial grid
fac = 1; % For Fourier interpolation
Nx = fac*Lx*4; Ny = fac*Ly*4;
x0 = 0; % Initial position of soli segment
y0 = 0;
dx = 2*Lx/Nx;
x = dx*[-Nx/2:Nx/2-1];
dy = 2*Ly/Ny;
y = dy*[-Ny/2:Ny/2-1];
[foo,xmin] = min(abs(x-(x0-220)));
[foo,xmax] = min(abs(x-(x0+120)));
[foo,ymin] = min(abs(y-(y0-200)));
[foo,ymax] = min(abs(y-(y0+200)));
xp = unique(round(linspace(xmin,xmax,250)));
yp = unique(round(linspace(ymin,ymax,250)));
yvec = linspace(y(ymin),y(ymax),400);

% Set up temporal grid
dt = 10; 
tout = [20,50]/dt;
Np = length(tout);


% Setup figure
f = 4;
fig_width = 13*4; % in cm
fig_height = 4.5*f;
ml = 1/13; % Margin left
mr = 0.02; % Margin right
mt = 0.04; % Margin top
mb = 0.22;  % Margin bottom
pb1 = 0.01; % Interaxes padding bottom
pr = 0.008; % Interaxes padding right
spanx1 = 3.3/13;
spanx2 = 1 - 2*spanx1-mr - 2*ml - pr;
spany = (1-mt-mb);

fontsize = 9*f;
cmap = load('CoolWarmFloat257.csv');

fh=figure(1);
clf();
fh.Renderer = 'Painters';
set(gcf,'Resize','off')
set(fh,'paperposition',[0,0,fig_width,fig_height],...
       'papersize',[fig_width,fig_height],'paperunits',...
       'centimeters','units','centimeters');

%Contour plots with dashed line exact solution
for ii=1:Np
    inc = tout(ii);
    load([data_dir,'/',sprintf('%05d',inc),'.mat'],'u');
%     u = u;
    % Interpolate if desired
    if fac > 1
        u = fftInterpolate(u,fac*size(u));
    end
    
    axes('Position',[ml+(ii-1)*(spanx1+pr),mb,spanx1,spany]);
    contourf(x(xp)-x0,y(yp),u(yp,xp),25,'edgecolor','none');
    colormap(cmap); %flipud(gray));%cmap);
    caxis(ran);
    set(gca,'FontSize',fontsize,'TickLabelInterpreter','latex',...
            'ytick',round([-200,0,200]),...
            'xtick',[-200,0,200]);%,'TickLength',[0.005 0.005]); 
    xlabel('$x$','interpreter','latex');
    if ii == 1
        ylabel('$y$','interpreter','latex');
    else
        set(gca,'yticklabel',{'','',''});
    end
    hold on
    [apred,qpred] = reg_exact_soln(yvec,inc*dt,a0,q0);
    xv = -cumtrapz(yvec,qpred);
    xtrue = xv-max(xv)+inc*dt*(a0/3+q0^2);
    plot(xtrue,yvec,':k','Linewidth',2.5*f)

    axis([x(xmin)-x0,round(x(xmax)-x0),round(y(ymin)),round(y(ymax))]);
end

set(fh,'paperposition',[0,0,fig_width,fig_height],...
       'papersize',[fig_width,fig_height],'paperunits',...
       'centimeters','units','centimeters');


   %% Amplitude Plot
  axes('Position',[2*ml+2*spanx1+pr,mb,spanx2,spany]);
  % Grid for analysis
xlow = -450;
ylow =  -500;
yhigh = 500;
[foo,xmin] = min(abs(x-(x0+xlow)));
[foo,ymin] = min(abs(y-(y0+ylow)));
[foo,ymax] = min(abs(y-(y0+yhigh)));
kx = pi/Lx*fftshift([-Nx/2:Nx/2-1]);
ky = pi/Ly*fftshift([-Ny/2:Ny/2-1]);
[KX,KY] = meshgrid(kx,ky);
fac = 1;

tout2 = [0 40 80 120 160];
Np = length(tout2); 

% Extract amplitude and slope as functions of y
 anumerics = zeros(Np,ymax-ymin+1);    

%Extract amplitude and angle from numerics
for ii=1:Np
    inc = tout2(ii)/dt;
    load([data_dir,'/',sprintf('%05d',inc),'.mat'],'u');
    
    % Interpolate if desired
    if fac > 1
        u = fftInterpolate(u,fac*size(u));
    end
 
  xmod = 2.3*inc*dt;
 [foo,xmax] = min(abs(x-xmod));
   
    v = u(ymin:ymax,xmin:xmax);
    [a1,I1] = max(v,[],2); 
   
    anumerics(ii,:)=a1;
end

amp = zeros(Np,ymax-ymin+1);
yvec = y(ymin:ymax);
 
%Exact solution
for jj = 1:length(tout2)
 [amp(jj,:),~]=reg_exact_soln(yvec,tout2(jj),a0,q0);
end 
   
hold on;
plot(yvec,ones(size(yvec)),'b-.',...
        'linewidth',f/2);
for jj=2:Np
plot(yvec,anumerics(jj,:),'k-',...
     'linewidth',f/2);
plot(yvec,amp(jj,:),'r--',...
        'linewidth',f/2);
end

set(gca,'FontSize',fontsize,'TickLabelInterpreter','latex',...
        'linewidth',f/2,'ytick',[0,1,2,3],'xtick',[-400,-200,0,200,400]);
xlabel('$y$','interpreter','latex');
ylabel('$a$','interpreter','latex');
axis([-500,500,0,2.4]);

hold off

text(25,.5,['$t = ',num2str(round(tout2(2))),'$'],'interpreter','latex',...
     'verticalalignment','middle','horizontalalignment','center',...
     'fontsize',7*f);
text(210,.5,['$',num2str(round(tout2(3))),'$'],'interpreter',...
     'latex',...
     'verticalalignment','middle','horizontalalignment','right',...
     'fontsize',7*f);
text(330,.5,['$',num2str(round(tout2(4))),'$'],'interpreter','latex',...
     'verticalalignment','middle','horizontalalignment','right',...
     'fontsize',7*f);
text(445,.5,['$',num2str(round(tout2(5))),'$'],'interpreter','latex',...
     'verticalalignment','middle','horizontalalignment','right',...
     'fontsize',7*f);

 doc_name = 'reg_ref_acc.pdf';  
  print(fh,'-dpdf',doc_name);
 %saveas(fh,['C:\Users\samry\Documents\KP\matlab_plots\mean_flow\',doc_name])


