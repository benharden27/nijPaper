clear

ranFL = 1; % 1 = whole year, % 2 = seasonal % 3= monthly;
savefold = '~/Documents/projects/kogur/NIJpaper/paper/figures/';
maxd = 800;
minx = -5; maxx = 69;

load ~/data/KGA/products/all_gridded.mat

% plotting setup
load '~/data/KGA/other/kogBathymetry'
x1 = 81.80;
regdist = regdist-x1;

tp = 5;

% set up
g = [1 1 1];
fontsz = 10;
load cm_kogurmovie
load cmaps_toBen.mat

 

PDvec = squeeze(nanmedian(PDfinal1(:,:,:),1));
Svec = squeeze(nanmedian(Sfinal1(:,:,:),1));
Tvec = squeeze(nanmedian(Tfinal1(:,:,:),1));
wsvec = squeeze(nanmean(ws(:,:,:),1));

[X,D] = meshgrid(xvec,dvec);
Svec = cfill(Svec);
Tvec = cfill(Tvec);
PDvec = cfill(PDvec);
wsvec = cfill2(wsvec);

% Tnew = Tvec;
% Tnew(PDvec<27.9) = nan;

sft = 0.03;

figure

% plot velocity
vran = -21:3:21;
mysub(3,1,1)
hold on
[C,h1] = contourf(xvec,dvec,wsvec*100,vran,'linestyle','none')
h2 = cbarf(wsvec*100,vran);
cmapmaker(h1,h2)
contour(xvec,dvec,wsvec,[0 0],'color',[.5 .5 .5])
% caxis([-.25 .25])
set(gca,'ydir','reverse','xdir','reverse')
% colormap(cm_kogurmovie)
axis([minx maxx 0 maxd])

PDcont = [27:.1:27.9 27.95 28 28 28.05];
[C,h] = contour(xvec,dvec,PDvec,PDcont,'k');
contour(xvec,dvec,PDvec,[27.8 27.8],'color','k','linewidth',2)
clabel(C,h,PDcont([1:8 10:length(PDcont)]))
% contour(xvec,dvec,Tnew,[0 0],'color','k','linewidth',2)

plot(xvec([1 end]),[650 650],'k--','linewidth',2)

ylabel(h2,'Current Velocity (cm s^-^1)')
ylabel('Depth (m)')
fill([regdist(~isnan(regdist))' regdist(end) regdist(1)],[regbat(~isnan(regdist))' 3000 3000],[.5 .5 .5])
plot(posxa(:,tp),posda(:,tp),'.','color',[.7 .7 .7],'markersize',10)

text(10,750,'Along Stream','horizontalalignment','center','backgroundcolor','w')

pos = get(gca,'position');
pos(3) = pos(3)-sft;
set(gca,'position',pos);
pos = get(h2,'position');
pos(1) = pos(1)-sft;
set(h2,'position',pos);

text(150,1000,'(a)','fontsize',14)

dxs = 2; dys = 20;
for i = 1:12
    text(dist(i),-40,num2str(i),'HorizontalAlignment','center')
    x = dist(i); y = 20;
    fill([x x-dxs x+dxs x],[y y-2*dys y-2*dys y],'k');
end
% for i = 1:12
%     text(dist(i),-25,num2str(i),'HorizontalAlignment','center')
% end
set(gca, 'Layer', 'top')

% Plot Temperature
mysub(3,1,2)
hold on
ccon = temp_levels(1:13)';
% ccon = [-2 -1 -.5 0 0.5 1 1.5 2 3 4];
% ccon = [-2 -1 -0.5 -0.25 0 0.25 0.5 0.75 1 1.5 2 3 4 5 6 8 10 12]';
[C,h1] = contourf(xvec,dvec,Tvec,ccon,'linestyle','none');
settempcm2(h1);
h2 = cbarf(Tvec,ccon','vert','nonlinear');
% cmapmaker(h1,h2,'ds_temp')
settempcm2(h2,Tvec);
ylabel(h2,'Potential Temperature (^oC)')
% set(h2,'fontsize',fontsz)

% PDcont = [22:.1:27.9 27.95 28 28.03 28.05];
[C,h] = contour(xvec,dvec,PDvec,PDcont,'k');
contour(xvec,dvec,PDvec,[27.8 27.8],'color','w','linewidth',2)
% contour(xvec,dvec,Tnew,[0 0],'color','k','linewidth',2)
plot(xvec([1 end]),[650 650],'k--','linewidth',2)
clabel(C,h,PDcont([1:8 10:length(PDcont)]))


set(gca,'ydir','reverse','xdir','reverse')
fill([regdist(~isnan(regdist))' regdist(end) regdist(1)],[regbat(~isnan(regdist))' 3000 3000],[.5 .5 .5])
plot(posx(:,tp),posd(:,tp),'.','color',[.7 .7 .7],'markersize',10)

axis([minx maxx 0 maxd])
ylabel('Depth (m)')
text(10,750,'Temperature','horizontalalignment','center','backgroundcolor','w')
dxs = 2; dys = 20;
for i = 1:12
    text(dist(i),-40,num2str(i),'HorizontalAlignment','center')
    x = dist(i); y = 20;
    fill([x x-dxs x+dxs x],[y y-2*dys y-2*dys y],'k');
end
pos = get(gca,'position');
pos(3) = pos(3)-sft;
set(gca,'position',pos);
pos = get(h2,'position');
pos(1) = pos(1)-sft;
set(h2,'position',pos);

text(150,1000,'(b)','fontsize',14)
set(gca, 'Layer', 'top')

% plot salinity

mysub(3,1,3)
hold on
ccon = sal_levels(2:13);
% ccon = [33 34 34.5 34.75 34.85 34.9 34.91 34.92 34.93];
[C,h1] = contourf(xvec,dvec,Svec,ccon,'linestyle','none');
h2 = cbarf(Svec,ccon,'vert','nonlinear');
% cmapmaker(h1,h2,'midy')

setsalcm2(h1);
setsalcm2(h2,Svec);
ylabel(h2,'Salinity')
% set(h2,'fontsize',fontsz)

% PDcont = [27:.1:27.9 27.95 28 28.03 28.05];
[C,h] = contour(xvec,dvec,PDvec,PDcont,'k');
contour(xvec,dvec,PDvec,[27.8 27.8],'color','w','linewidth',2)
% contour(xvec,dvec,Tnew,[0 0],'color','k','linewidth',2)
clabel(C,h,PDcont([1:8 10:length(PDcont)]))
% plot(xvec([1 end]),[650 650],'color','k','linewidth',2)

plot(xvec([1 end]),[650 650],'k--','linewidth',2)



set(gca,'ydir','reverse','xdir','reverse')
fill([regdist(~isnan(regdist))' regdist(end) regdist(1)],[regbat(~isnan(regdist))' 3000 3000],[.5 .5 .5])
plot(posx(:,tp),posd(:,tp),'.','color',[.7 .7 .7],'markersize',10)

axis([minx maxx 0 maxd])
xlabel('Distance (km)')
ylabel('Depth (m)')
text(10,750,'Salinity','horizontalalignment','center','backgroundcolor','w')
set(gcf,'renderer','painters')

pos = get(gca,'position');
pos(3) = pos(3)-sft;
set(gca,'position',pos);
pos = get(h2,'position');
pos(1) = pos(1)-sft;
set(h2,'position',pos);

text(150,1000,'(c)','fontsize',14)
dxs = 2; dys = 20;
for i = 1:12
    text(dist(i),-40,num2str(i),'HorizontalAlignment','center')
    x = dist(i); y = 20;
    fill([x x-dxs x+dxs x],[y y-2*dys y-2*dys y],'k');
end
% set(gcf,'paperunits','normalized','paperposition',[0 0 1 1])
% print('-depsc',[savefold 'mean_section'])
set(gca, 'Layer', 'top')

print_fig('plot_section',savefold,0.82,1);

%     cropall(savefold);


