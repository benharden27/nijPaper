clear
load ~/Documents/projects/kogur/NIJpaper/code/wavespecs tvec cv fs wp wn
savefold = '~/Documents/projects/kogur/NIJpaper/paper/figures/';
load cm_midwhite 
load invgray
% cv = wsmean(:,3)*100;
% gamma=3;beta=2;
% fs=morsespace(gamma,beta,{0.05,pi},pi/1000,4);
% [wp,wn]=wavetrans(cv,conj(cv),{gamma,beta,fs,'bandpass'});
% h=wavespecplot(tvec,cv,2./fs,sqrt(squared(wp)+squared(wn)));

yran = 0:200;

glb = 0.1;
gap = 0.05;
sw = 0.15;
mw = 1-2*glb-sw-gap;
mh = (1-glb-2*gap)/2;


figure
subplot('position',[glb mh+gap+glb mw mh]), hold on
plot(tvec,real(cv),'k')
plot(tvec,imag(cv),'color',[.5 .5 .5])
set(gca,'xtick',datenum(2011,9:21,1),'xlim',tvec([1 end]),'xticklabel','','box','on','tickdir','out')
% datetick('x',19,'keepticks')
% set(gca,)
ylabel('Velocity (cm/s)')

subplot('position',[glb glb mw mh]), hold on
contourf(tvec,2./fs,sqrt(squared(wp)+squared(wn))',0:35,'linestyle','none')
ylabel('Period (days)')
set(gca,'xtick',datenum(2011,9:21,1))
datetick('x',12,'keepticks')
set(gca,'xlim',tvec([1 end]),'ylim',[1.5 300],'ydir','reverse','yscale','log','tickdir','out','box','on')
plot(tvec([1 end]),[8 8],'k--')
% cbarf(sqrt(squared(wp)+squared(wn)),0:35)
pos = get(gca,'position');
h = colorbar;
set(h,'position',[glb+mw+gap glb+mh+gap 0.03 mh])
ylabel(h,'Amplitude (cm/s)')
set(gca,'position',pos)

subplot('position',[glb+mw+gap glb sw mh]), hold on
plot(mean(sqrt(squared(wp)+squared(wn))),2./fs,'k')
set(gca,'ylim',[1.5 300],'ydir','reverse','yaxislocation','right','yscale','log','tickdir','out','box','on')
plot(get(gca,'xlim'),[8 8],'k--')
ylabel('Period (days)')
xlabel('Mean Amplitude (cm/s)')
% packfig(2,5,'rows')

hAllAxes = findobj(gcf,'type','axes');


colormap(invgray)
% pos = get(gca,'position');
% hc = colorbar;
% hcpos = get(hc,'position');
% hcpos(1) = hcpos(1) + 0.075;
% set(hc,'position',hcpos)
% ylabel(hc,'Magnitude (cm/s)')
% set(gca,'position',pos)

print_fig('TRWwavelet',savefold,1,0.5)
