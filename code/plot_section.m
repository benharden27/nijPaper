clear

% Select some user parameters
ranFL = 1; % 1 = whole year, % 2 = seasonal % 3= monthly;
savefold = '~/Documents/projects/kogur/NIJpaper/paper/figures/';
maxd = 800;
minx = -5; maxx = 69;

% load in the gridded data
load ~/data/KGA/products/all_gridded.mat

% Load in bathymetry
load '~/data/KGA/other/kogBathymetry'
x1 = 81.80;
regdist = regdist-x1;

% which timestep has the locations of all data well represented?
tp = 5;

% set up other plotting parameters
fontsz = 10;
% load cmaps_toBen.mat
temp_levels = [ -2.0000    -1.0000   -0.5000   -0.2500  0 ...
                0.2500    0.5000    0.7500    1.0000    1.5000 ...
                2.0000    3.0000    4.0000    5.0000    6.0000 ...
                8.0000   10.0000   12.0000];
sal_levels = [31.0000   33.0000    34.0000   34.2000   34.4000...
              34.6000   34.8000   34.8500   34.8800   34.9000...
              34.9100   34.9200   34.9400   34.9700   35.0000...
              35.1000   35.2000   35.3000];

            
% process data to produce means/medians
PDvec = squeeze(nanmedian(PDfinal1(:,:,:),1));
Svec = squeeze(nanmedian(Sfinal1(:,:,:),1));
Tvec = squeeze(nanmedian(Tfinal1(:,:,:),1));
wsvec = squeeze(nanmean(ws(:,:,:),1));

% fill data to edge of domain
[X,D] = meshgrid(xvec,dvec);
Svec = cfill(Svec);
Tvec = cfill(Tvec);
PDvec = cfill(PDvec);
wsvec = cfill2(wsvec);

% The shift in page dimensions to accommodate cbar
sft = 0.03;

% begin figure
figure

% plot velocity
mysub(3,1,1)
hold on

% set velocity range
vran = -14:2:14;

% fill contour the velocity field
[C1,h1] = contourf(xvec,dvec,wsvec*100,vran,'linestyle','none');
h2 = cbarf(vran,vran); % add colorbar
contour(xvec,dvec,wsvec,[0 0],'color',[.5 .5 .5]) % add zero line
set(gca,'ydir','reverse','xdir','reverse')
axis([minx maxx 0 maxd])

% overlay the density field
PDcont = [27:.1:27.9 27.95 28 28 28.05];
[C,hPD] = contour(xvec,dvec,PDvec,PDcont,'k');
contour(xvec,dvec,PDvec,[27.8 27.8],'color','k','linewidth',2)
clabel(C,hPD,PDcont([1:8 10:length(PDcont)]))

% add depth of sill
plot(xvec([1 end]),[650 650],'k--','linewidth',2)

% label axes
ylabel(h2,'Current Velocity (cm s^-^1)')
ylabel('Depth (m)')

% add bathymetry
fill([regdist(~isnan(regdist))' regdist(end) regdist(1)],[regbat(~isnan(regdist))' 3000 3000],[.5 .5 .5])

% plot instrument locations
plot(posxa(:,tp),posda(:,tp),'.','color',[.7 .7 .7],'markersize',10)

% add figure labels
text(10,650,'Along Stream Velocity','horizontalalignment','center','backgroundcolor','w')
text(150,1000,'(a)','fontsize',14)

% move colorbar and axes for good spacing
pos = get(gca,'position');
pos(3) = pos(3)-sft;
set(gca,'position',pos);
pos = get(h2,'position');
pos(1) = pos(1)-sft;
set(h2,'position',pos);

% add mooring locations and labels
dxs = 2; dys = 20;
for i = 1:12
    text(dist(i),-40,['KGA ' num2str(i)],'HorizontalAlignment','center')
    x = dist(i); y = 20;
    fill([x x-dxs x+dxs x],[y y-2*dys y-2*dys y],'k');
end

% Move axes to top
set(gca, 'Layer', 'top')

% change colormaps for velocity
pause(.01) % For some reason, matlab needs a break at this point
[~,cm2] = cmapmaker(C1,h1,h2);



% Plot Temperature
mysub(3,1,2)
hold on
ccon = temp_levels(1:13)';
[~,h1] = contourf(xvec,dvec,Tvec,ccon,'linestyle','none');
h2 = cbarf(Tvec,ccon','vert','nonlinear');
setTScm(h2,'temp');
ylabel(h2,'Potential Temperature (^oC)')

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
text(10,650,'Temperature','horizontalalignment','center','backgroundcolor','w')
dxs = 2; dys = 20;
for i = 1:7
    text(dist(i),-40,['KGA ' num2str(i)],'HorizontalAlignment','center')
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

pause(0.01)
setTScm(h1,'temp');


% plot salinity

mysub(3,1,3)
hold on
ccon = sal_levels(5:13);
% ccon = [33 34 34.5 34.75 34.85 34.9 34.91 34.92 34.93];
[C,h1] = contourf(xvec,dvec,Svec,ccon,'linestyle','none');
h2 = cbarf(ccon,ccon,'vert','nonlinear');
% cmapmaker(h1,h2,'midy')

setTScm(h2,'sal')
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
text(10,650,'Salinity','horizontalalignment','center','backgroundcolor','w')
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
    text(dist(i),-40,['KGA' num2str(i)],'HorizontalAlignment','center')
    x = dist(i); y = 20;
    fill([x x-dxs x+dxs x],[y y-2*dys y-2*dys y],'k');
end
% set(gcf,'paperunits','normalized','paperposition',[0 0 1 1])
% print('-depsc',[savefold 'mean_section'])
set(gca, 'Layer', 'top')

pause(0.01)
setTScm(h1,'sal');


print_fig('plot_section',savefold,0.82,1);

%     cropall(savefold);


