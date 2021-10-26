%% Plots line soliton and Y soliton 

xl = -128;
xr = 128;
yl = -128;
yr = 128;

xvec = linspace(xl,xr,1.5*(xr-xl));
yvec = linspace(yl,yr,1.5*(yr-yl));

[X,Y] = meshgrid(xvec,yvec);

Np = 2;

% Setup figure
ml = 0.15; % Margin left
mr = 0.1; % Margin right
mt = 0.04; % Margin top
mb = 0.22;  % Margin bottom
pb = 0.01; % Interaxes padding bottom
pr = 0.07; % Interaxes padding right
spanx = (1-ml-mr-(Np-1)*pr)/Np;
spany = (1-mt-mb);
f = 4; % Factor to increase figure size (dashed line hack)
fig_width = 12*f; % in cm
fig_height = 4.5*f;
fontsize = 9*f;
cmap = load('CoolWarmFloat257.csv');

fh=figure(1);
clf();
fh.Renderer = 'Painters';
set(gcf,'Resize','off')
set(fh,'paperposition',[0,0,fig_width,fig_height],...
       'papersize',[fig_width,fig_height],'paperunits',...
       'centimeters','units','centimeters');

   
%Plot of line soliton
axes('Position',[ml,mb,spanx,spany]);
u = sech((X+.9*Y)./sqrt(12)).^2;
hold on
contourf(X,Y,u,25,'edgecolor','none');
colormap(cmap); %flipud(gray));%cmap);
caxis([0,1]);
plot(zeros(size(yvec)),yvec,'w--','linewidth',f/2)
    circular_arrow(fh,89,[0 0], 110,36,-1,'w',f*2)
text(-33,105,'$\varphi$','interpreter',...
      'latex',...
      'verticalalignment','middle','horizontalalignment','center',...
      'fontsize',10*f,'color','white');
text(95,105,'$(a)$','interpreter','latex',...
    'verticalalignment','middle','horizontalalignment','left',...
      'fontsize',10*f,'color','white')
hold off
    
set(gca,'FontSize',fontsize,'TickLabelInterpreter','latex',...
         'ytick',[],...
            'xtick',[]);%,'TickLength',[0.005 0.005]); 
    xlabel('$x$','interpreter','latex');
    ylabel('$y$','interpreter','latex');
    axis([xl,xr,yl,yr]);
set(fh,'paperposition',[0,0,fig_width,fig_height],...
       'papersize',[fig_width,fig_height],'paperunits',...
       'centimeters','units','centimeters');
   
%Plot of Y soliton
axes('Position',[ml+spanx+pr,mb,spanx,spany]);
k1 = -1;
k2 = 0;
k3 = 1.5;
a1 = (k2-k1)^2/4;
a2 = (k3-k1)^2/4;
a3 = (k3-k2)^2/4;
uf = @(z1,z2) (4*a1^2*a2.*z1+4*a2*a3^2.*z2+...
            (4*a1^2*a3+8*sqrt(a1^3*a3^3)+4*a1*a3^2).*z1.*z2)...
            ./(a2+a1*z1+a3*z2).^2;
z1_i =  exp(sqrt(1/12)*((k1-k2).*X+(k1^2-k2^2).*Y/2));
z2_i =  exp(sqrt(1/12)*((k3-k2).*X+(k3^2-k2^2).*Y/2));
uy = uf(z1_i,z2_i);

contourf(X,Y,uy,25,'edgecolor','none')
axis([xl,xr,yl,yr]);
text(-10,90,'$(a_2,q_2)$','interpreter',...
      'latex',...
      'verticalalignment','middle','horizontalalignment','left',...
      'fontsize',10*f,'color','white');
text(53,-48,'$(a_3,q_3)$','interpreter',...
      'latex',...
      'verticalalignment','middle','horizontalalignment','left',...
      'fontsize',10*f,'color','white');
text(-50,-48,'$(a_1,q_1)$','interpreter',...
      'latex',...
      'verticalalignment','middle','horizontalalignment','right',...
      'fontsize',10*f,'color','white');
text(-125,105,'$(b)$','interpreter','latex',...
    'verticalalignment','middle','horizontalalignment','left',...
      'fontsize',10*f,'color','white')
xlabel('$x$','interpreter','latex');

set(gca,'FontSize',fontsize,'TickLabelInterpreter','latex',...
         'ytick',[],...
            'xtick',[]);
        
 doc_name = 'soliton_solns.pdf';  
  print(fh,'-dpdf',doc_name);        

