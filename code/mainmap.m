clear
savefold = '~/Documents/projects/kogur/NIJpaper/paper/figures/';

fontsz = 14;
fold = 'mat_36lp_20m';
dir = '~/data/KGA/ADCP/';
dir = [dir fold '/'];
files = dirr([dir '/*.mat']);
nfile = size(files,1);
% nfile = 7;


dir = '~/data/KGA/ADCP/final/';
for i = 1:nfile
    if i == 11
        files = dirr([dir 'KGA' num2str(i+1) '_*.mat']);
    else
        files = dirr([dir 'KGA' num2str(i) '_*.mat']);
    end
    file = files(1,:);
    if file(end) == ' '
        file = file(1:end-1);
    end
    load([dir file])

    poslat(i) = str2num(data.meta.latitude(1:end-2));
    poslon(i) = -str2num(data.meta.longitude(1:end-2));
   
end

distm = sw_dist(poslat,poslon,'km');

%%



load ~/data/bathymetry/bathy_ibcao_ds

lon1 = findnear(lon,-32);
lon2 = findnear(lon,-15);
lat1 = findnear(lat,65);
lat2 = findnear(lat,79);

lambda = mean(lat([lat1 lat2]));
dx = sw_dist([lambda lambda],lon([lon1 lon2]));
dy = sw_dist(lat([lat1 lat2]),lon([lon1 lon1]));
ratio = dx/dy;

xscale = sw_dist([lambda lambda],[1 2]);
yscale = sw_dist([1 2],lon([lon1 lon1]));

st = 5;
load invgray


figure, hold on
colormap(invgray)

contourf(lon(lon1:lon2),lat(lat1:lat2),-bath(lat1:lat2,lon1:lon2),[0:250:2000],'linestyle','none')
contour(lon(lon1:lon2),lat(lat1:lat2),-bath(lat1:lat2,lon1:lon2),[0 0],'k')

caxis([0 3000])
h = cbarf([0 nanmax(nanmax(-bath(lat1:lat2,lon1:lon2)))],[0:250:2000]);
ylabel(h,'Depth (m)')
pos_main = get(gca,'position');
pos_cb = get(h,'position');
set(h,'position',[pos_cb(1)-.05 pos_main(2) pos_cb(3) pos_main(4)])
set(gca,'plotboxaspectratio',[ratio 1 1],'box','on')
hold on


% load '/Users/benjamin/data/ERAinterim/N128/coastline'
% plot(coastline(:,1),coastline(:,2),'k','linewidth',1)


% plot currents
dotFL = 0;
lwd = 8;
NIJcol = [2,56,88]/255;
EGCcol = [44 162 95]/255;
NIICcol = [239,101,72]/255;

xpred = [-18.5:-0.1:-25.5];
xbrk = -21 ;
xi = findnear(xpred,xbrk);

NIJx = [-18.2081 -18.9923 -20.6318 -22.2179 -23.2872 -23.6970 -24.2317 -24.7307 -25.2118];
NIJy = [68.2705 67.8465 67.7764 67.7296 67.5360 67.4025 67.1555 66.9218 66.7282];
NIJs = spline(NIJx,NIJy);
ypred = ppval(NIJs,xpred);
% ypred(1:5:xi) = nan;
% ypred(2:5:xi) = nan;

plot(xpred(1:xi),ypred(1:xi),'color',NIJcol,'linestyle','--','linewidth',lwd)
plot(xpred(xi:end),ypred(xi:end),'color',NIJcol,'linewidth',lwd)
% plot(xpred,ypred,NIJcol,'linewidth',lwd)
if(dotFL)
    plot(NIJx,NIJy,'k.','markersize',10);
end
sc = 20;
tht = 45;
fill(xpred(end-1) - [sc/2*cosd(90-tht) sc*cosd(tht) -sc/2*cosd(90-tht)]/xscale,...
    ypred(end-1) - [-sc/2*sind(90-tht) sc*sind(tht) sc/2*sind(90-tht)]/yscale,...
    NIJcol,'linestyle','none')


ypred = 66.8:0.05:71.5;
EGCx = [-19.3102 -18.7005 -18.9819 -20.2482 -21.7489 -24.2876 -25.2193 -25.7821 -26.1104];
EGCy = [71.1008 70.2575 69.5548 69.1683 68.9575 68.5007 68.0615 67.3939 66.8668];
EGCs = spline(EGCy,EGCx);
plot(ppval(EGCs,ypred),ypred,'color',EGCcol,'linewidth',lwd)
if(dotFL)
plot(EGCx,EGCy,'k.','markersize',10);
end
xpred = ppval(EGCs,ypred);
tht = 75;
fill(xpred(2) - [sc/2*cosd(90-tht) sc*cosd(tht) -sc/2*cosd(90-tht)]/xscale,...
    ypred(2) - [-sc/2*sind(90-tht) sc*sind(tht) sc/2*sind(90-tht)]/yscale,...
    EGCcol,'linestyle','none')




ypred = 66.95:0.05:71.5;
ybrk = 68.25;
yi = findnear(ypred,ybrk);

EGC2x = [-19.3102 -18.7005 -18.6819 -19.8261 -21.2330 -22.8745 -24.0000 -24.7035 -25.4069];
EGC2y = [71.1008 70.2575 69.5548 68.8521 68.3777 67.9034 67.5520 67.2182 66.8317];
EGC2s = spline(EGC2y,EGC2x);
xpred = ppval(EGC2s,ypred);
plot(xpred(1:yi),ypred(1:yi),'color',EGCcol,'linewidth',lwd)
plot(xpred(yi:end),ypred(yi:end),'color',EGCcol,'linestyle','--','linewidth',lwd)
if(dotFL)
plot(EGC2x,EGC2y,'k.','markersize',10);
end
tht = 60;
fill(xpred(2) - [sc/2*cosd(90-tht) sc*cosd(tht) -sc/2*cosd(90-tht)]/xscale,...
    ypred(2) - [-sc/2*sind(90-tht) sc*sind(tht) sc/2*sind(90-tht)]/yscale,...
    EGCcol,'linestyle','none')


xpred = [-28:0.05:-20];
NIICx = [-27.2989 -26.4468 -24.6487 -23.3789 -21.7777 -19.8729];
NIICy = [64.7838 65.6215 66.5626 67.1210 67.4882 67.5105];
NIICs = spline(NIICx,NIICy);
ypred = ppval(NIICs,xpred);
% ypred(1:5:xi) = nan;
% ypred(2:5:xi) = nan;

plot(xpred,ypred,'color',NIICcol,'linewidth',lwd)
if(dotFL)
    plot(NIICx,NIICy,'k.','markersize',10);
end
tht = 170;
fill(xpred(end-1) - [sc/2*cosd(90-tht) sc*cosd(tht) -sc/2*cosd(90-tht)]/xscale,...
    ypred(end-1) - [-sc/2*sind(90-tht) sc*sind(tht) sc/2*sind(90-tht)]/yscale,...
    NIICcol,'linestyle','none')

fs = 14;
text(-20.5,70.3,'EGC','color',EGCcol,'fontsize',fs,'fontweight','bold')
text(-25,69,'sbEGC','color',EGCcol,'fontsize',fs,'fontweight','bold')
text(-20.2,68.45,'sEGC','color',EGCcol,'fontsize',fs,'fontweight','bold')
text(-21.5,68,'NIJ','color',NIJcol,'fontsize',fs,'fontweight','bold')
text(-21.5,67.2,'NIIC','color',NIICcol,'fontsize',fs,'fontweight','bold')
text(-28,66.3,{'Denmark','Strait','Sill'},'color','k','fontsize',fs-2,'HorizontalAlignment','center')
text(-28,69.3,'Greenland','color','k','fontsize',fs-2,'HorizontalAlignment','center')
text(-19,65.3,'Iceland','color','k','fontsize',fs-2,'HorizontalAlignment','center')

plot(poslon,poslat,'k.','markersize',10)
plot(poslon(1:7),poslat(1:7),'k.','markersize',20)

plot(-18-50/60,68,'k+','markersize',10,'linewidth',3)



print_fig('mainmap',savefold,1,.6)
