clear

savefold = '~/Documents/projects/kogur/NIJpaper/paper/figures/';

minlon = -30;
maxlon = -18;
minlat = 65;
maxlat = 70;

load ~/Documents/projects/kogur/NIJpaper/paper/code/TRWwaveparms

% calculate splines
% load ~/data/WHOI/bathymetry/bathy
% bath = elev; lon = lons; lat = lats';

load ~/data/bathymetry/bathy_ibcao_ds

bath = bath(find(lat>minlat,1)-1:find(lat>maxlat,1),find(lon>minlon,1)-1:find(lon>maxlon,1));
lon = lon(find(lon>minlon,1)-1:find(lon>maxlon,1));
lat = lat(find(lat>minlat,1)-1:find(lat>maxlat,1));
lon1 = lon; lat1 = lat;

[m,n] = size(bath);

for i = 1:n
    Y(:,i) = my_lowpass(bath(:,i),60);
end

for i = 1:m
    Y(i,:) = my_lowpass(Y(i,:),60);
end

figure
load invgray
colormap(invgray)
[c,h] = contourf(lon,lat,-Y,[0:250:2000],'k','linestyle','none');
caxis([0 3000])
h = cbarf([1 max(-Y(:))],[0:250:1750]);
ylabel(h,'Depth (m)')
pos_main = get(gca,'position');;
pos_cb = get(h,'position');
set(h,'position',[pos_cb(1) pos_main(2) pos_cb(3) pos_main(4)])

hold on
contour(lon,lat,bath,[0 0],'k')



% load ~/data/KGA/products/all_gridded.mat lon lat

prd = 4;

for i = 2:5
    loc = [poslon(i) poslat(i)];
    k = kvec(i,:);
%     k = nanmean(kvec(2:4,:));
    [CG,WCHECK,LOC,K,BSLOPE,LAMBDA,H] = trwtrace_nij(loc,k,prd,parms);
    LOC(H<-1.3*2,:) = nan;
    plot(LOC(:,1),LOC(:,2),'k')
    H2(i,:,:) = H;
    LOC2(i,:,:) = LOC;
    CG2(i,:,:) = CG;
    K2(i,:,:) = K;
    BSLOPE2(i,:,:) = BSLOPE;
    LAMBDA2(i,:,:) = LAMBDA;
end


plot(poslon,poslat,'k.','markersize',20)

axis([-26 -22 66.7 68.2])

dx = sw_dist([mean(get(gca,'ylim')) mean(get(gca,'ylim'))],get(gca,'xlim'));
dy = sw_dist(get(gca,'ylim'),poslon([1 1]));
ratio = dx/dy;
xscale = sw_dist([mean(get(gca,'ylim')) mean(get(gca,'ylim'))],[1 2]);
yscale = sw_dist([1 2],poslon([1 1]));

set(gca,'plotboxaspectratio',[ratio 1 1])




% clear LK6 KK6 LOCK6
% for i = 2:5
%     ii = find(H2(i,:)<-2,1);
%     LK6(i) = 64*LAMBDA2(i,ii)/LAMBDA2(i,1);
%     KK6(i,1) = kvec(i,1)*K2(i,ii,1)/K2(i,1,1);
%     KK6(i,2) = kvec(i,2)*K2(i,ii,2)/K2(i,1,2);
%     LOCK6(i,:) = LOC2(i,ii,:);
%     
% %     plot(LOCK6(i,1),LOCK6(i,2),'k.','markersize',20)
%     quiver(LOCK6(i,1),LOCK6(i,2),100000*KK6(i,1)/xscale,100000*KK6(i,2)/yscale,'color','k')
%     quiver(poslon(i),poslat(i),100000*kvec(i,1)/xscale,100000*kvec(i,2)/yscale,'color','k')
% 
% end



for i = 1:7
%    labm{i} = ['KGA' num2str(i)]; 
   labm{i} = ['KGA-' num2str(i)]; 
   dlon(i) = 0.05;
   switch(i)
       case{1}
           dlon(i) =  0.05;
           dlat(i) =  -0.02;
       case{2}
           dlat(i) = 0.02;
           dlon(i) = 0.07;
       otherwise
           dlat(i) = 0.03;
   end
end

text(poslon+dlon,poslat+dlat,labm)


print_fig('wave_trace',savefold,1,.6)