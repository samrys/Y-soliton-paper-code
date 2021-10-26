%% Creates schematic of Mach reflection and reg reflection

% Setup figure
ml = 0.09; % Margin left
mr = 0.02; % Margin right
mt = 0.04; % Margin top
mb = 0.21;  % Margin bottom
pb = 0.01; % Interaxes padding bottom
pr = 0.04; % Interaxes padding right
Nplots = 2;
spanx = (1-ml-mr-(Nplots-1)*pr)/Nplots;
spany = (1-mt-mb-(1-1)*pb)/1;
f = 4; % Factor to increase figure size (dashed line hack)
fig_width = 12*f; % in cm
fig_height = 4.5*f;
fontsize = 9*f;


 fh=figure(1);
 clf();
 fh.Renderer = 'Painters';
 set(fh,'Resize','off')
 set(fh,'paperposition',[0,0,fig_width,fig_height],...
       'papersize',[fig_width,fig_height],'paperunits',...
       'centimeters','units','centimeters');
% 

tl = linspace(-2,0,300);
tr = linspace(0,3.72,300);

%% Left plot: 
  axes('Position',[ml,mb,spanx,spany]);
  hold on;
%   axis([-3,5,-4,4]);
  
  %Represents the bottom
  plot(tl,zeros(size(tl)),'k-','linewidth',f)
  plot(tr,zeros(size(tr)),'k--','linewidth',f/4)
  plot(tr,.25*tr,'k-','linewidth',f)

  %Represents the first line soliton
  soli1 = linspace(0,2.8,150);
  plot(-1*ones(size(soli1)),soli1,'k-','linewidth',f/2)
  
  %Represents the second line soliton
  soli2 = linspace(.6,2.8,150);
  plot(1*ones(size(soli2)),soli2,'k-','linewidth',f/2)
  y1 = linspace(.25,.6,100);
  x1 = -.15*(y1-.6)+1;
  plot(x1,y1,'k-','linewidth',f/1.25)
  yy1 = linspace(.6,1.06,100);
  xx1 = .8+(yy1-.9)*-2/3-.05*(yy1-.9).*(yy1-.6);
  plot(xx1,yy1,'k-','linewidth',f/2);

  
    %Represents the third line soliton
  soli3 = linspace(1.5,2.8,150);
  plot(3*ones(size(soli3)),soli3,'k-','linewidth',f/2)
  y2 = linspace(.75,1.5,100);
  x2 = -.15*(y2-1.5)+3;
  plot(x2,y2,'k-','linewidth',f/1.25)
  yy2 = linspace(1.5,2.65,100);
  xx2 = 2.8+(yy2-1.9)*-.5-.1*(yy2-1.9).*(yy2-1.5);
  plot(xx2,yy2,'k-','linewidth',f/2);

  
  %Bottom boundary
%   xz = linspace(-3,-2.9,50);
%   yz = -xz-3;
%   plot(xz,yz,'k-','linewidth',f/4)
%   
 annotation('arrow',[.226 .264],[.77 .77],'linewidth',f/2);
 annotation('arrow',[.316 .356],[.77 .77],'linewidth',f/2);
 circular_arrow(fh,2.5,[0 0],6,12,-1,'k',f);
 text(2.85,.41,'$\varphi$','interpreter',...
      'latex',...
      'verticalalignment','middle','horizontalalignment','center',...
      'fontsize',7*f);
 hold off;  
 set(gca,'FontSize',fontsize,'TickLabelInterpreter','latex',...
      'XTick',[],'YTick',[],'YTickLabel',[],'XTickLabel',[],'visible','off');
 axis([-3.5,5.5,-4,4]);
 text(-2,3.5,'$(a)$','interpreter','latex',...
    'verticalalignment','top','horizontalalignment','left',...
      'fontsize',10*f,'color','black')
 set(fh,'paperposition',[0,0,fig_width,fig_height],...
       'papersize',[fig_width,fig_height],'paperunits',...
       'centimeters','units','centimeters');

   
   
%% Right plot
 axes('Position',[ml+spanx+pr,mb,spanx,spany]);
  hold on;
%   axis([-3,5,-4,4]);
  tl = linspace(-2.8,0,300);
  tr = linspace(0,2.8,300);
  %Represents the bottom
  plot(tl,zeros(size(tl)),'k-','linewidth',f/4)
  plot(tr,zeros(size(tr)),'k-','linewidth',f/4)
%   plot(tr,-.5*tr,'k-','linewidth',f)

  %Represents the first line soliton
  soli1 = linspace(0,2.8,150);
  plot(.28*soli1-1.5,soli1,'k-','linewidth',f/2)
  plot(.28*soli1-1.5,-soli1,'k-','linewidth',f/2)
  
  %Dashed line for reflection
  plot(-1.5*ones(size(soli1)),soli1,'k--','linewidth',f/4)
  plot(-1.5*ones(size(soli1)),-soli1,'k--','linewidth',f/4)
  
  %Represents the second line soliton
  a = 0.1;
  b = .32;
  fd = -3;
  c = -0.1;
  d = 1;
  x1b = linspace(a,c,100);
  y1b = quadrat(a,b,fd,c,d,x1b);
  soli2 = linspace(.32,2.8,150);
  plot(.28*(soli2-soli2(1)+.25),soli2,'k-','linewidth',f/2)
  plot(.28*(soli2-soli2(1)+.25),-soli2,'k-','linewidth',f/2)
  y1 = linspace(-.34,.34,200);
  x1 = 0.09*ones(size(y1));
  
%   y1b = linspace(.25,.587,100);
%   x1b = .253+(y1b-.587)*(.33-.253)/(.25-.587)+-.6*(y1b-.25).*(y1b-.587);
  plot(x1b,y1b,'k-','linewidth',f/2)
  plot(x1b,-y1b,'k-','linewidth',f/2)
  plot(x1,y1,'k-','linewidth',f/1.25)
  
    %Represents the third line soliton
  a = 2.05;
  b = .8;
  fd = -3;
  c = 1.4;
  d = 2.2;
  x3 = linspace(a,c,100);
  y3 = quadrat(a,b,fd,c,d,x3);
  soli3 = linspace(.8,2.8,150);
  plot(.28*soli3+1.826,soli3,'k-','linewidth',f/2)
  plot(.28*soli3+1.826,-soli3,'k-','linewidth',f/2)
  y2 = linspace(-.8,.8,100);
  x2 = 2.05*ones(size(y2));
  plot(x2,y2,'k-','linewidth',f/1.25)
  plot(x3,y3,'k-','linewidth',f/2)
  plot(x3,-y3,'k-','linewidth',f/2)
  
  %Bottom boundary
%   xz = linspace(-3,-2.9,50);
%   yz = -xz-3;
%   plot(xz,yz,'k-','linewidth',f/4)
%   
 annotation('arrow',[.71 .74],[.725 .7],'linewidth',f/2);
 annotation('arrow',[.645 .675],[.795 .77],'linewidth',f/2);
 circular_arrow(fh,2.1,[-1.5 0],82,14,1,'k',f);
 circular_arrow(fh,2.1,[-1.5 0],278,14,-1,'k',f);
 text(-1.15,2.5,'$\varphi$','interpreter',...
      'latex',...
      'verticalalignment','middle','horizontalalignment','center',...
      'fontsize',7*f);
 text(-1.2,-2.5,'$-\varphi$','interpreter',...
      'latex',...
      'verticalalignment','middle','horizontalalignment','center',...
      'fontsize',7*f); 
 hold off;  
 set(gca,'FontSize',fontsize,'TickLabelInterpreter','latex',...
      'XTick',[],'YTick',[],'YTickLabel',[],'XTickLabel',[],'visible','off');
 axis([-2.5,6.5,-4,4]);
 set(fh,'paperposition',[0,0,fig_width,fig_height],...
       'papersize',[fig_width,fig_height],'paperunits',...
       'centimeters','units','centimeters');

text(-2.5,3.5,'$(b)$','interpreter','latex',...
    'verticalalignment','top','horizontalalignment','left',...
      'fontsize',10*f,'color','black')
% 
print(fh,'-dpdf','mach_reflection_prsa.pdf');


function [ q ] = quadrat(a,b,fd,c,d,x)
    K1 = (d-b)./(c-a);
    K2 = (fd - K1)./(a-c);
    q =  b + (x-a).*K1 + (x-a).*(x-c).*K2;
 end     
