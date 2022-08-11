function [] = display()
set(0,'defaulttextinterpreter','latex')                                   ;%
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

T     = readtable('./dat/parameters.csv','PreserveVariableNames',true);
par   = table2array(T)
nx    = par(1)
ny    = par(2)
nsave = par(7)
T     = readtable('./dat/x.csv','PreserveVariableNames',true);
D     = table2array(T)
x     = D(:,1)
T     = readtable('./dat/y.csv','PreserveVariableNames',true);
D     = table2array(T)
y     = D(:,1)
T     = readtable('./dat/zhs.csv','PreserveVariableNames',true);
D     = table2array(T)
hs    = D(:,2)
hs    = reshape(hs,nx,ny);

load('./colormaps/devon');
hmax = 1e-1
v = VideoWriter(['runoff' '.mp4'],'MPEG-4');v.FrameRate=30;v.Quality=100;open(v);
hold on
for sim=1:nsave
    % load data
    T = readtable(['./dat/tdt_',num2str(sim),'.csv'],'PreserveVariableNames',true);
    D = table2array(T);
    t = D(1);
    T = readtable(['./dat/hQxQy_',num2str(sim-1),'.csv'],'PreserveVariableNames',true);
    D = table2array(T);
    h = reshape(sqrt(D(:,2).^2+D(:,3).^2),nx,ny);
    h = reshape(D(:,1),nx,ny);
    h(h<1e-2)=0.0;
    
fig1=figure(1);
clf;
%plot first data 
ax1 = axes; 
im = imagesc(hs'); 
im.AlphaData = 0.75; % change this value to change the background image transparency 
colormap(ax1,gray) 
    xlabel('Easting [m]');
    ylabel('Northing [m]');
    axis equal;
axis tight;
t = seconds(t);
t.Format = 'hh:mm:ss';

title(['$t=',char(t),'$'],'FontSize',14);
set(gca,'FontSize',12,'TickLabelInterpreter','latex');
hold all; 
%plot second data 
ax2 = axes; 
im1 = imagesc(h');
alpha = max(h'./hmax,hmax);
im1.AlphaData = alpha; % change this value to change the foreground image transparency 
colormap(ax2,flip(devon)) 
caxis(ax2,[0 hmax])
set(gca,'FontSize',12,'TickLabelInterpreter','latex');
axis equal;
axis tight;
%link axes 
linkaxes([ax1,ax2]) 
ax2.Visible = 'off'; 
ax2.XTick = []; 
ax2.YTick = [];
set(ax2,'color','none','visible','off');
%add differenct colormap to different data if you wish 
%set the axes and colorbar position  
cb2 = colorbar(ax2,'FontSize',10,'TickLabelInterpreter','latex','Fontsize',10,'Location','southoutside','Position',[0.3750 0.10 0.25 0.0625]);  
        cb2.Label.String     ='$h(x,y)|_{t}$ [m]';
        cb2.Label.Units      ='normalized';
        %cb.Label.Position   =[0.5 2.2];
        cb2.Label.FontSize   =10;
        cb2.Label.Interpreter='Latex';
set(gcf,'color','white');
drawnow

        F = getframe(gcf);
        writeVideo(v,F.cdata);
end
close(v);

end