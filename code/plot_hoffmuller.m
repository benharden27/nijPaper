clear

savefold = '~/Documents/projects/kogur/NIJpaper/paper/figures/';
load ~/data/KGA/products/all_gridded.mat

%% set up parameters

maxx = 70;
maxxi = find(xvec>maxx,1);
xvec = xvec(1:maxxi);
xloc = unique(posx(~isnan(posx)));

PD = PDfinal1(:,:,1:maxxi);
for i = 1:length(xvec)
    for j = 1:length(tvec)
        PDin = squeeze(PD(j,:,i));
        goodi = ~isnan(PDin);
        d28(j,i) = interp1(PDin(goodi),dvec(goodi),27.8);
    end
end


g = 9.8;
Pvec = sw_pres(dvec,67);
for xi = 1:length(xvec)
    S = Sfinal1(:,:,xi);
    T = Tfinal1(:,:,xi);
    P = repmat(Pvec,[length(tvec) 1]);
    T = sw_temp(S,T,P,zeros(size(P)));
    GA(:,:,xi) = sw_gpan(S',T',P');
end
GA = permute(GA,[2 1 3]);
D = (GA(findnear(dvec,300),:)/g);

for i = 1:length(xvec)
    gpan(:,:,i) = sw_gpan(Sfinal1(:,:,i),Tfinal1(:,:,i),P);
end



umean = squeeze(nanmean(ugrid(:,:,1:maxxi),2));
vmean = squeeze(nanmean(vgrid(:,:,1:maxxi),2));
wsmean = squeeze(nanmean(ws(:,:,1:maxxi),2));
PDmean = squeeze(PDfinal1(:,5,1:maxxi));
GPAN = squeeze(GA(:,5,:));

for i = 1:length(xvec)
    lsti = find(~isnan(GA(100,:,i)),1,'last');
    GPAN(:,i) = GA(:,lsti,i);
end

for i = 1:length(xvec)
    uh(:,i) = my_bandpass(umean(:,i),3*2,3*8);
    vh(:,i) = my_bandpass(vmean(:,i),3*2,3*8);
    wsh(:,i) = my_bandpass(wsmean(:,i),3*2,3*8);
    wsl(:,i) = my_lowpass(wsmean(:,i),3*8);
    d28h(:,i) = my_bandpass(d28(:,i),3*2,3*8);
    PDh(:,i) = my_bandpass(PDmean(:,i),3*2,3*8);
    GPANh(:,i) = my_bandpass(GPAN(:,i),3*2,3*8);
end
d28h(:,1) = nan;

EKE = sqrt(uh.^2+vh.^2);
for i = 1:length(xvec)
    EKE(:,i) = my_lowpass(EKE(:,i),3*8);
end

load ~/Documents/projects/kogur/NIJpaper/code/wavespecs_amp POW

for i = 1:length(xvec)
    [th,x,y] = varelip(uh(:,i),vh(:,i),0);
    ang(i)=-th;
    wsr(:,i) = rot_ac(uh(:,i),vh(:,i),-th+90+180);
end


%% plot correlations


%% subset of parms
% n = 3;
% gap = 0.02;
% mgap = gap*4;
% bot1 = 0.05;
% lef = 0.05;
% wid = 0.9;
% hei = (1-(gap*5+mgap+bot1))/(n*2);
% for i = 1:n*2
%     bot(i) = bot1 + (i-1)*(hei+gap);
% end
% bot(n+1:end) = bot(n+1:end) + mgap;
% bot = bot(end:-1:1);
% 
% ti = findnear(tvec,datenum(2012,3,1));
% 
% figure
% for i = 1:3
%     for j = 1:2
%         subplot('position',[lef bot((j-1)*n+i) wid hei]), hold on
%         switch i
%             case 1
%                 varplot = wsmean'*100;
%                 levs = -20:2:20;
%             case 2
%                 varplot = wsr'*100;
%                 levs = -20:2:20;
%             case 3
% %                 varplot = EKE';
% %                 levs = 0:0.005:0.2;
%                 varplot = POW;
%                 levs = [0:1:20];
%         end
%     [~,h1] = contourf(tvec,xvec,varplot,levs,'linestyle','none');
%     cmapmaker(h1);
%     set(gca,'ydir','reverse','xtick',datenum(2011,9:21,1))
%     datetick('x',19,'keepticks')
%     switch(j)
%         case 1
%             axis([tvec([1 ti]) xvec([1 end])])
%         case 2
%             axis([tvec([ti end]) xvec([1 end])])
%     end
%             
%     switch i
%         case 1
%             title('Barotropic along stream mean')
%         case 2
%             title('8 day high-passed (rotated to var elips)')
%         case 3
%             title('Wavelet amplitude at 4 days')
%             tplot = xvec/39;
%             tplots = datenum(2011,9,1:10:1000);
%             for k = 1:length(tplots)
%                 plot(tplots(k)+tplot,xvec(end:-1:1),'k')
%             end
%     end
%     end
% %     hlines(xloc)
% end
% 
% print_fig('hoff_vars',savefold,1,0.7)
% 


%% subset of parms
n = 3;
gap = 0.02;
mgap = gap*4;
bot1 = 0.05;
lef = 0.05;
wid = 0.8;
hei = (1-(gap*5+mgap+bot1))/(n*2);
for i = 1:n*2
    bot(i) = bot1 + (i-1)*(hei+gap);
end
bot(n+1:end) = bot(n+1:end) + mgap;
bot = bot(end:-1:1);

ti = findnear(tvec,datenum(2012,3,1));

tline = datenum(2011,9:21,1);

figure
for i = 1:3
        switch i
            case 1
                subplot(1,10,1:3), hold on
                varplot = wsmean'*100;
                
                levs = -20:2:20;
            case 2
                subplot(1,10,4:6), hold on
                varplot = wsr'*100;
                levs = -20:2:20;
            case 3
%                 varplot = EKE';
%                 levs = 0:0.005:0.2;
                subplot(1,10,7:9), hold on
                varplot = POW;
                levs = [-20:2:20];
        end
        varplot(varplot<-20) = -20;
                varplot(varplot>20) = 20;
    [~,h1] = contourf(xvec,tvec,varplot',levs,'linestyle','none');
    cmapmaker(h1);
    set(gca,'ytick',datenum(2011,9:21,1))
    datetick('y','mmmyy','keepticks')
    axis([xvec([1 end]) tvec([1 end])])
    xlabel('Distance (km)')
    hlines(tline,'k--')
    if(i>1)
        set(gca,'yticklabel',[])
    end
%     pos = get(gca,'position');
%     set(gca,'position',[pos(1) pos(2)+0.1 pos(3) pos(4)-0.1])
    switch i
        case 1
            title('(a)')
            pos = get(gca,'position');
            h2 = cbarf(varplot,levs);
            cmapmaker(h1,h2)
            set(h2,'position',[0.87 pos(2) 0.03 pos(4)])
            ylabel(h2,'cm s^-^1')
            set(gca,'position',pos)
        case 2
            title('(b)')
            
        case 3
            title('(c)')
            tplot = xvec/20;
            tplots = datenum(2011,9,1:10:1000);
            for k = 1:length(tplots)
                plot(xvec(end:-1:1),tplots(k)+tplot,'k')
            end
            
    end
    
    
%     hlines(xloc)
end

print_fig('hoff_vars',savefold,1,.9)




%% Energy correlation
% figure
% % col = {'k','k','k','k',[.75 .75 .75]};
% % sty = {'-','--','-.',':','-'};
% for i = 1:5
%     [varcorr,lags] = myxcorr(POW(i+1,:),POW(i,:),0,0,1/8);
% %     varcorr = varcorr(end:-1:1);
% %     lags = lags(end:-1:1);
% %     midi = (length(varcorr)+1)/2;
% %     maxlag(i) = lags(midi - 1 + find(varcorr(midi:end)==max(varcorr(midi:end))));
%     plot(lags,varcorr,'k')
% %     phalag=360*lags/4;
%     lagoff(i) = lags(find(varcorr==max(varcorr)));
%     hold on
%     
%     
% end
% if i == 5
%     legend('KGA2-KGA3','KGA3-KGA4','KGA4-KGA5','KGA5-KGA6','KGA6-KGA7')
% else
%     legend('1-2','2-3','3-4','4-5')
% end 
% axis([-360 360 -1 1])
% plot(get(gca,'xlim'),[0 0],'k')
% plot([0 0],[-1 1],'k')
% 



