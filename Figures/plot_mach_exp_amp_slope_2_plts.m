% Extract amplitude and slope from bent soliton evolution
data_dir = 'C:\Users\samry\Documents\MATLAB\KP2\data\miles_resonant\v_shape_a_1_q_p7';
load([data_dir,'\parameters.mat'],'Lx','Ly');

% Set up spatial grid
fac = 4; % For Fourier interpolation
Nx = fac*Lx*4; Ny = fac*Ly*4;
a0 = 1;
q0 = .7;
x0 = 0; % Initial position of center
y0 = 0;
dx = 2*Lx/Nx;
x = dx*[-Nx/2:Nx/2-1];
dy = 2*Ly/Ny;
y = dy*[-Ny/2:Ny/2-1];
% Grid for analysis
xlow = -360;
xhigh = 600;
ylow =  -500;
yhigh = 500;
[foo,xmin] = min(abs(x-(x0+xlow)));
[foo,xmax] = min(abs(x-(x0+xhigh)));
[foo,ymin] = min(abs(y-(y0+ylow)));
[foo,ymax] = min(abs(y-(y0+yhigh)));
kx = pi/Lx*fftshift([-Nx/2:Nx/2-1]);
ky = pi/Ly*fftshift([-Ny/2:Ny/2-1]);
[KX,KY] = meshgrid(kx,ky);

xvec = linspace(xlow,xhigh,10000);



% Set up temporal grid for plotting
dt = 10; 
tout = [0 50 150 250];
Np = length(tout);

% Setup figure
ml = 0.09; % Margin left
mr = 0.02; % Margin right
mt = 0.03; % Margin top
mb = 0.17;  % Margin bottom
pb = 0.068; % Interaxes padding bottom
pr = 0.085; % Interaxes padding right
Nplots = 2;
top = .22;
bot = 1-mt-mb-pb-top; 
spanx = (1-ml-mr-(Nplots-1)*pr)/Nplots;
%spany = (1-mt-mb-(1-1)*pb)/1;
f = 4; % Factor to increase figure size (dashed line hack)
fig_width = 12*f; % in cm
fig_height = 5.5*f;
fontsize = 9*f;

% Extract amplitude and slope as functions of y
 anumerics1 = zeros(Np,ymax-ymin+1);    %Stem and right-most legs
 qnumerics1 = zeros(Np,ymax-ymin+1);            
 anumerics2 = anumerics1;               %Stem and left-most legs
 qnumerics2 = qnumerics1;

%Extract amplitude and angle from numerics
for ii=1:Np
    inc = tout(ii)/dt;
    load([data_dir,'/',sprintf('%05d',inc),'.mat'],'u');
    
    % Interpolate if desired
    if fac > 1
        u = fftInterpolate(u,fac*size(u));
    end
    
    v = u(ymin:ymax,xmin:xmax);
    m1 = islocalmax(v,2,'MaxNumExtrema',2,'MinProminence',.02); %Finds local maxima of each row
    v2 = v;
    v2(m1~=1)=NaN;
    [a1,I1] = max(v2,[],2); %Right-most legs
    [a2,I2] = min(v2,[],2); %Left-most legs
    
    a2(I1==I2) = NaN; %Eliminate duplicates
    x1 = x(xmin+I1-1);
    x2 = x(xmin+I2-1);
    x2(I1==I2) = NaN;
    
    anumerics1(ii,:)=a1;
    anumerics2(ii,:)=a2;
    qnumerics1(ii,:)=x1;
    qnumerics2(ii,:)=x2;
end

%Smooth for q and take derivative of y with respect to x
for ii=1:Np
    qnumerics1(ii,:) = smooth(qnumerics1(ii,:),fac*10);
    qnumerics2(ii,:) = smooth(qnumerics2(ii,:),fac*40);
end
qnumerics1(1:end,2:end-1) = -(qnumerics1(1:end,3:end)-...
                             qnumerics1(1:end,1:end-2))/(2*dx);
                        
qnumerics2(1:end,2:end-1) = -(qnumerics2(1:end,3:end)-...
                             qnumerics2(1:end,1:end-2))/(2*dx);

for ii=1:Np
    qnumerics1(ii,:) = smooth(qnumerics1(ii,:),fac*3);
    qnumerics2(ii,:) = smooth(qnumerics2(ii,:),fac*2);
end


%Remove values that don't make sense
 qnumerics1(qnumerics1>.8)=nan;
 qnumerics1(qnumerics1<-.8)=nan;
 
 qnumerics2(qnumerics2<-1.7)=nan;
 qnumerics2(qnumerics2>1.7)=nan;
  
 qnumerics2(isnan(anumerics2))=NaN;

%Initialize for the exact solution
amp1 = zeros(Np,ymax-ymin+1);
amp2 = amp1;
ang1 = amp1;
ang2 = amp1;
yvec = y(ymin:ymax);
 
%Exact solution
for jj = 1:length(tout)
 [amp1(jj,:),amp2(jj,:),ang1(jj,:),ang2(jj,:)]=mach_exact_soln(yvec,tout(jj),a0,q0);
end

a_init = a0.*ones(size(yvec));
q_init = (q0.*(yvec<=1)-q0.*(yvec>0)).*ones(size(yvec));
an_init = a0.*ones(size(yvec));
qn_init = q0*tanh(-yvec./1);

% 
% xshift = -20;
% xmat = repmat(xvec,5,1);
% xmat = xmat-(amp<.43).*xshift;

fh=figure(1);
clf();
fh.Renderer = 'Painters';
set(fh,'paperposition',[0,0,fig_width,fig_height],...
       'papersize',[fig_width,fig_height],'paperunits',...
       'centimeters','units','centimeters');


   %% Creating plot
   % Amplitude
clf
axes('Position',[ml,mb,spanx,bot]);
hold on;
plot(yvec,ones(size(yvec))*(sqrt(a0)+q0)^2,'b-.',...
        'linewidth',.75*f);
for jj=2:Np
plot(yvec,anumerics1(jj,:),'k-',...
     'linewidth',f/2);
plot(yvec,amp1(jj,:),'r--',...
        'linewidth',f/2);
end


plot(yvec,an_init,'k-','linewidth',f/2);
plot(yvec,a_init,'r--','linewidth',f/2);
hold off;

set(gca,'FontSize',fontsize,'TickLabelInterpreter','latex',...
        'linewidth',f/2,'ytick',[0,1,2,3],'xtick',[-50,0,50]);
xlabel('$y$','interpreter','latex');
ylabel('$a$','interpreter','latex');
axis([-72,72,0,3.2]);
text(0,.97,['$t = ',num2str(round(tout(1))),'$'],'interpreter','latex',...
     'verticalalignment','top','horizontalalignment','center',...
     'fontsize',7*f);
text(1,1.7,['$t = ',num2str(round(tout(2))),'$'],'interpreter',...
     'latex',...
     'verticalalignment','middle','horizontalalignment','center',...
     'fontsize',7*f);
text(25,1.7,['$',num2str(round(tout(3))),'$'],'interpreter','latex',...
     'verticalalignment','middle','horizontalalignment','center',...
     'fontsize',7*f);
text(50,1.7,['$',num2str(round(tout(4))),'$'],'interpreter','latex',...
     'verticalalignment','middle','horizontalalignment','right',...
     'fontsize',7*f);


%Simple wave inset
axes('Position',[ml,mb+pb+bot,spanx,top]);
hold on;
plot(yvec,q0^2.*ones(size(yvec)),'b-.','linewidth',.75*f);
for jj=2:Np
plot(yvec,anumerics2(jj,:),'k-',...
     'linewidth',f/2);
plot(yvec,amp2(jj,:),'r--',...
        'linewidth',f/2);
end



set(gca,'FontSize',fontsize,'TickLabelInterpreter','latex',...
        'linewidth',f/2,'ytick',[0,.5,1],'xtick',[200,400]);

ylabel('$a_{\rm p}$','interpreter','latex');
axis([0,460,0,.7]);
text(15,.25,['$t = ',num2str(round(tout(2))),'$'],'interpreter',...
     'latex',...
     'verticalalignment','middle','horizontalalignment','left',...
     'fontsize',7*f);
text(230,.25,['$t=',num2str(round(tout(3))),'$'],'interpreter','latex',...
     'verticalalignment','middle','horizontalalignment','center',...
     'fontsize',7*f);
text(435,.25,['$t=',num2str(round(tout(4))),'$'],'interpreter','latex',...
     'verticalalignment','middle','horizontalalignment','right',...
     'fontsize',7*f);
hold off

% Angle
axes('Position',[ml+spanx+pr,mb,spanx,bot]);
hold on;
plot(yvec,zeros(size(yvec)),'b-.','linewidth',.75*f);
for jj=2:Np
plot(yvec,qnumerics1(jj,:),'k-',...
     'linewidth',f/2);
plot(yvec,ang1(jj,:),'r--',...
        'linewidth',f/2);
end
plot(yvec,qn_init,'k-','linewidth',f/2);
plot(yvec,q_init,'r--','linewidth',f/2);


set(gca,'FontSize',fontsize,'TickLabelInterpreter','latex',...
        'linewidth',f/2,'ytick',[-1,0,1],...
        'xtick',[-50,0,50]);
xlabel('$y$','interpreter','latex');
ylabel('$q$','interpreter','latex');
axis([-72,72,-1,1.2]);
text(.05,.6,['$t=',num2str(round(tout(1))),'$'],'interpreter','latex',...
     'verticalalignment','middle','horizontalalignment','left',...
     'fontsize',7*f);
text(-12,.2,['$',num2str(round(tout(2))),'$'],'interpreter',...
     'latex',...
     'verticalalignment','middle','horizontalalignment','right',...
     'fontsize',7*f);
text(-45,.2,['$',num2str(round(tout(3))),'$'],'interpreter','latex',...
     'verticalalignment','middle','horizontalalignment','left',...
     'fontsize',7*f);
text(-77,.2,['$t = ',num2str(round(tout(4))),'$'],'interpreter','latex',...
     'verticalalignment','middle','horizontalalignment','left',...
     'fontsize',7*f); 
 
hold off;
 %Simple wave inset
axes('Position',[ml+spanx+pr,mb+pb+bot,spanx,top]);
hold on;
plot(yvec,sqrt(a0)*ones(size(yvec)),'b-.','linewidth',.75*f);
for jj=2:Np
plot(yvec,qnumerics2(jj,:),'k-',...
     'linewidth',f/2);
plot(yvec,ang2(jj,:),'r--',...
        'linewidth',f/2);
end


set(gca,'FontSize',fontsize,'TickLabelInterpreter','latex',...
        'linewidth',f/2,'ytick',[0,1,2],'xtick',[200,400]);
ylabel('$q_{\rm p}$','interpreter','latex');
axis([0,460,.5,2]);

text(20,1.3,['$t = ',num2str(round(tout(2))),'$'],'interpreter',...
     'latex',...
     'verticalalignment','middle','horizontalalignment','left',...
     'fontsize',7*f);
text(220,1.3,['$t =',num2str(round(tout(3))),'$'],'interpreter','latex',...
     'verticalalignment','middle','horizontalalignment','center',...
     'fontsize',7*f);
text(447,1.3,['$ t = ',num2str(round(tout(4))),'$'],'interpreter','latex',...
     'verticalalignment','middle','horizontalalignment','right',...
     'fontsize',7*f);

% 
%   print(fh,'-dpdf','mach_amp_slope_a_1_q_p7_2plts.pdf');
