% Plot results from mean flow soliton simulation
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
[foo,xmin] = min(abs(x-(x0-180)));
[foo,xmax] = min(abs(x-(x0+420)));
[foo,ymin] = min(abs(y-(y0-245)));
[foo,ymax] = min(abs(y-(y0+245)));
xp = unique(round(linspace(xmin,xmax,250)));
yp = unique(round(linspace(ymin,ymax,250)));

% Set up temporal grid
dt = 10; 
% t = linspace(1,24,25)*dt;
tout = [0,100,250]/dt;
Np = length(tout);
% tout_inds = zeros(1,Np);
% for ii=1:Np
%     [foo,ind] = min(abs(t-tout(ii)));
%     tout_inds(ii) = ind;
% end
% tout = t(tout_inds);


% Setup figure
ml = 0.11; % Margin left
mr = 0.03; % Margin right
mt = 0.04; % Margin top
mb = 0.22;  % Margin bottom
pb = 0.01; % Interaxes padding bottom
pr = 0.01; % Interaxes padding right
spanx = (1-ml-mr-(Np-1)*pr)/Np;
spany = (1-mt-mb);
fig_width = 12; % in cm
fig_height = 4.5;
fontsize = 9;
cmap = load('CoolWarmFloat257.csv');


fh=figure(1);
clf();
fh.Renderer = 'Painters';
set(gcf,'Resize','off')
set(fh,'paperposition',[0,0,fig_width,fig_height],...
       'papersize',[fig_width,fig_height],'paperunits',...
       'centimeters','units','centimeters');

for ii=1:Np
    inc = tout(ii);
    load([data_dir,'/',sprintf('%05d',inc),'.mat'],'u');
%     u = u;
    % Interpolate if desired
    if fac > 1
        u = fftInterpolate(u,fac*size(u));
    end
    
    axes('Position',[ml+(ii-1)*(spanx+pr),mb,spanx,spany]);
    contourf(x(xp)-x0,y(yp),u(yp,xp),25,'edgecolor','none');
    colormap(cmap); %flipud(gray));%cmap);
    caxis(ran);
    set(gca,'FontSize',fontsize,'TickLabelInterpreter','latex',...
            'ytick',round([-200,0,200]),...
            'xtick',[-200,0,200,400]);%,'TickLength',[0.005 0.005]); 
    xlabel('$x$','interpreter','latex');
    if ii == 1
        ylabel('$y$','interpreter','latex');
    else
        set(gca,'yticklabel',{'','',''});
    end
    axis([x(xmin)-x0,round(x(xmax)-x0),round(y(ymin)),round(y(ymax))]);
    
    %Texts
    if ii==1
       text(165,200,'$(a_0,q_0)$','interpreter',...
      'latex',...
      'verticalalignment','top','horizontalalignment','left',...
      'fontsize',fontsize,'color','white');
       text(165,-200,'$(a_0,-q_0)$','interpreter',...
      'latex',...
      'verticalalignment','bottom','horizontalalignment','left',...
      'fontsize',fontsize,'color','white');
    elseif ii==3
       text(220,0,'$(a_{\rm i},q_{\rm i})$','interpreter',...
      'latex',...
      'verticalalignment','middle','horizontalalignment','right',...
      'fontsize',fontsize,'color','white');       
       text(65,200,'$(a_{\rm n},q_{\rm n})$','interpreter',...
      'latex',...
      'verticalalignment','top','horizontalalignment','right',...
      'fontsize',fontsize,'color','white');       
       text(65,-200,'$(a_{\rm n},-q_{\rm n})$','interpreter',...
      'latex',...
      'verticalalignment','bottom','horizontalalignment','right',...
      'fontsize',fontsize,'color','white'); 
    end
end


set(fh,'paperposition',[0,0,fig_width,fig_height],...
       'papersize',[fig_width,fig_height],'paperunits',...
       'centimeters','units','centimeters');

 doc_name = ['v_shape_a_1_q_p7_',num2str(tout(1)*dt),...
     '_',num2str(tout(2)*dt),'_',num2str(tout(3)*dt),'.pdf'];  
 print(fh,'-dpdf',doc_name);
 %saveas(fh,['C:\Users\samry\Documents\KP\matlab_plots\mean_flow\',doc_name])


