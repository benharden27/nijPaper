%% Load data from moorings

clear

savefold = '~/Documents/projects/kogur/NIJpaper/paper/figures/';

fontsz = 14;
fold = 'mat_36lp_20m';
dir = '~/data/KGA/ADCP/';
dir = [dir fold '/'];
files = dirr([dir '/*.mat']);
nfile = size(files,1);
nfile = 7;

% load ~/data/KGA/products/all_gridded.mat

for i = 1:nfile
    load([dir files(i,:)])
%     nadcp(i) = str2num(files(i,4));
    nadcp(i,:) = files(i,4:5);
    eval(['u' nadcp(i,:) '= u;'])
    eval(['v' nadcp(i,:) '= v;'])
    eval(['d' nadcp(i,:) '= d;'])
    eval(['t' nadcp(i,:) '= t;'])
    clear u v t d
end

vars = char(who('t*'));
nlast = str2num(vars(nfile,2:3));
varstr = [vars repmat(' ',[size(vars,1) 1])]';
varstr = varstr(:)';
dvec = eval(['d' vars(nfile,2:3)]);
dvec = 15:5:900;
dd = diff(dvec([1 2]));
dt = 1/24;
tvec = eval(['min([' varstr ']):dt:max([' varstr '])']);
xvec = [1:nfile];
% xvec = [1 3:5];
unew = nan(nfile,length(tvec),length(dvec));
vnew = nan(nfile,length(tvec),length(dvec));

for i = 1:nfile
    % add on depths
    di = eval(['findnear(dvec,d' nadcp(i,:) '(1))']);
    df = eval(['findnear(dvec,d' nadcp(i,:) '(end))']);
    ti = eval(['findnear(tvec,t' nadcp(i,:) '(1))']);
    tf = eval(['findnear(tvec,t' nadcp(i,:) '(end))']);
%     if i==8|i==9;
%         tf = tf-1;
%     end
    unew(i,ti:tf,di:df) = eval(['u' nadcp(i,:) '''']);
    vnew(i,ti:tf,di:df) = eval(['v' nadcp(i,:) '''']);
end
   

dioffb = 38; %distance off bottom
dioffb = 20; %distance off bottom

for i = 1:nfile
    dbot = find(~isnan(unew(i,3000,:)),1,'last')-dioffb;
%     dbot(dbot>38) = 38;
%     dbot(dbot>18) = 18;
%     dbot = 18;
    ubot1(i,:) = squeeze(unew(i,:,dbot));
    vbot1(i,:) = squeeze(vnew(i,:,dbot));
end

ubot2 = ubot1(:,51:7000); vbot2 = vbot1(:,51:7000);

ubot = ubot1(:,51:7000);
vbot = vbot1(:,51:7000);
tbot = tvec(51:7000);

lowcut = 2;
highcut = 8;
for i = 1:nfile
    ubot(i,:) = my_bandpass(ubot(i,:),lowcut*24,highcut*24);
    vbot(i,:) = my_bandpass(vbot(i,:),lowcut*24,highcut*24);
end

di = findnear(dvec,100);
umean = squeeze(nanmean(unew(:,:,di:end),3));
vmean = squeeze(nanmean(vnew(:,:,di:end),3));
% wsmean = squeeze(nanmean(ws,2));
for i = 1:length(xvec)
    uh(:,i) = my_bandpass(umean(i,:),lowcut*24,highcut*24);
    vh(:,i) = my_bandpass(vmean(i,:),lowcut*24,highcut*24);
%     wsh(:,i) = my_bandpass(wsmean(:,i),lowcut*3,highcut*3);
end

% uh = ubot';
% vh = vbot';

clear wsr
for di = 1:nfile
    [th,x,y] = varelip(uh(:,di),vh(:,di),0);
    % principle axis angle (angle is counter-clockwise from east)
    ang(di)=-th;
    wsr(:,di) = rot_ac(uh(:,di),vh(:,di),-th+90+180);
end

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


%% set up bathymetric requirements

load ~/data/bathymetry/bathy_ibcao_ds
dlon = diff(lon([1 2]));
dlat = diff(lat([1 2]));
glo = 1.5;
gloi = round(glo/dlon);
gla = .3;
glai = round(gla/dlon);

lon1 = findnear(lon,poslon(end))-gloi;
lon2 = findnear(lon,poslon(1))+gloi;
lat1 = findnear(lat,poslat(1))-glai;
lat2 = findnear(lat,poslat(end))+glai;

lambda = mean(lat([lat1 lat2]));
dx = sw_dist([lambda lambda],lon([lon1 lon2]));
dy = sw_dist(lat([lat1 lat2]),lon([lon1 lon1]));
ratio = dx/dy;

xscale = sw_dist([lambda lambda],[1 2]);
yscale = sw_dist([1 2],lon([lon1 lon1]));

[m,n] = size(bath);
for i = 1:n
    Y(:,i) = my_lowpass(bath(:,i),5);
end
for i = 1:m
    Y(i,:) = my_lowpass(Y(i,:),5);
end

%% plot the phase correlation



figure
col = {'k','k','k','k',[.75 .75 .75]};
sty = {'-','--','-.',':','-'};
for i = 1:5
    [varcorr,lags] = myxcorr(wsr(:,i+1),wsr(:,i+2),0,0,1);
    varcorr = varcorr(end:-1:1);
    lags = lags(end:-1:1);
    midi = (length(varcorr)+1)/2;
    maxlag(i) = lags(midi - 1 + find(varcorr(midi:end)==max(varcorr(midi:end))));
    plot(360*lags/4,varcorr,'color',col{i},'linestyle',sty{i},'linewidth',2)
    phalag=360*lags/4;
    phaoffi(i) = phalag(find(varcorr==max(varcorr)));
    hold on
    
    
end
if i == 5
    legend('KGA2-KGA3','KGA3-KGA4','KGA4-KGA5','KGA5-KGA6','KGA6-KGA7')
else
    legend('1-2','2-3','3-4','4-5')
end 
axis([-360 360 -1 1])
plot(get(gca,'xlim'),[0 0],'k')
plot([0 0],[-1 1],'k')



%% calculate the TRW predicitons

dx = sw_dist([mean(poslat([1 5])) mean(poslat([1 5]))],poslon([1 5]),'km');
dy = sw_dist(poslat([1 5]),poslon([1 1]),'km');
arrayang = -atand(dx/dy);

T = 3.6*24; % measured time period of waves (hours)
Tsec = T*60*60;
phaoff = mean(phaoffi(1:3));
Del = mean(ang(2:4))-arrayang; % angle between phase propegation and mooring line (deg)
Ds = mean(distm(2:4)); % average instrument spacing (km)

% Phase velocity = 

cp = (1/T)*(360/phaoff)*(Ds/cosd(Del)); % in km/hour

for i = 1:5000
    rand1 = randi(3);
    rand2 = randi(3);
    rand3 = randi(3);
    cpmc(i) = (1/T)*(360/phaoffi(rand1))*(Ds/cosd(ang(rand3+1)-arrayang))*24;
end

cpms = cp * 1e3 / 3600;
cpkd = cp * 24;
dcp = cpkd*sqrt((1.25/phaoff)^2+(0.2/Ds)^2+(0.03/1)^2);

lam = cp*T;
dlam = lam*dcp/cpkd;
% 

% Calculate vertical stability
load ~/data/KGA/products/all_gridded.mat PDfinal1 dvec xvec
PD = squeeze(nanmean(PDfinal1));
P = repmat(sw_pres(dvec,68),[length(xvec) 1])';
[~,gradrhoz]  = gradient(PD,xvec,dvec);
bvf2 = 9.8.*gradrhoz./PD;
N = nanmean(bvf2(2:end,:));
N = nanmean(N(1:7));


% Calculate coriolis frequency
OM = 7.2921e-5;
lambda = mean(poslat(2:4));
f = 2*OM*sind(lambda);

% Bathymetric angle and bottom slope
D = 500;
xll = [0 cumsum(sw_dist(repmat(poslat(3),[length(lon) 1]),lon,'km'))'];
yll = [0 cumsum(sw_dist(lat,repmat(poslon(3),[length(lat) 1]),'km'))'];
[Zx,Zy] = gradient(bath,xll,yll);
Zxm = mean(mean(Zx(findnear(lat,poslat(2)):findnear(lat,poslat(4)),findnear(lon,poslon(4)):findnear(lon,poslon(2)))));
Zym = mean(mean(Zy(findnear(lat,poslat(2)):findnear(lat,poslat(4)),findnear(lon,poslon(4)):findnear(lon,poslon(2)))));
angB = atand(Zxm/Zym);
Ta = sqrt(Zxm.^2+Zym.^2)/1000;      

% angle between phase propagation and downslope
th = Del+arrayang-angB;

% adjuct lam and create k
lam = lam*1e3;
k = 2*pi/lam;

% calculate theorotical period
Tcalc = 2*pi*tanh(2*pi*N*D/lam/f)/N/Ta/sind(th)/24/60/60;
Tsec2 = Tcalc*24*60*60;

% Calculate the theoretic angle to downslope
thcalc = asind(2*pi*tanh(2*pi*N*D/lam/f)/Tsec/N/Ta);


Cgx = N*Ta*cosd(th)^2/k/tanh(N*D*k/f) - N^2*Ta*D*sind(th)^2/f/sinh(N*D*k/f)^2;
Cgy = -N*Ta*sind(2*th)/2/k/tanh(N*D*k/f) - N^2*Ta*D*sind(2*th)/2/f/sinh(N*D*k/f)^2;
Cg = sqrt(Cgx^2+Cgy^2);
Cgkd = Cg*3600*24/1e3;
Cgang = atand(-Cgy/Cgx)+90;
Dtrap = f/N/k;


display('Results of the TRW calculation')
display('--------------')
display(['Wave period = ' num2str(T/24)])
display(['Range of phase offset, KGA2-5 = ' num2str(phaoffi(1:3))])
display(['Average phase offset, KGA2-5 = ' num2str(mean(phaoffi(1:3)))])
display(['Range of angle to array, KGA2-5 = ' num2str(ang(2:4)-arrayang)])
display(['Average angle to array, KGA2-5 = ' num2str(Del)])
display(['Phase speed = ' num2str(cpkd)])
display(['Phase speed error = ' num2str(dcp)])
display(['Wavelength = ' num2str(lam/1000)])
display(['Wavelength error = ' num2str(dlam)])
display('--------------')
display(['Average BVF = ' num2str(N)])
display(['Coriolis frequency = ' num2str(f)])
display(['Bottom Slope = ' num2str(Ta)])
display(['Calculated angle from downslope = ' num2str(thcalc)])
display(['Measured angle from downslope = ' num2str(th)])
display('--------------')
display(['Theoretical Group Speed = ' num2str(Cgkd)])
display(['Theoretical Group Propegation = ' num2str(Cgang)])
display(['Trapping Depth = ' num2str(Dtrap)])

clear kvec
for i = 1:length(poslon)
    kvec(i,1) = k*sind(ang(i));
    kvec(i,2) = k*cosd(ang(i));
end

clear parms
parms(1) = f;
parms(2) = N;
parms(3) = D;
save ~/Documents/projects/kogur/NIJpaper/code/TRWwaveparms poslat poslon kvec parms

%% plot the positions
% load the bathymetry and scale
load ~/data/bathymetry/bathy_ibcao_ds
dlon = diff(lon([1 2]));
dlat = diff(lat([1 2]));
glo = 1.5;
gloi = round(glo/dlon);
gla = .3;
glai = round(gla/dlon);

lon1 = findnear(lon,poslon(end))-gloi;
lon2 = findnear(lon,poslon(1))+gloi;
lat1 = findnear(lat,poslat(1))-glai;
lat2 = findnear(lat,poslat(end))+glai;

lambda = mean(lat([lat1 lat2]));
dx = sw_dist([lambda lambda],lon([lon1 lon2]));
dy = sw_dist(lat([lat1 lat2]),lon([lon1 lon1]));
ratio = dx/dy;

xscale = sw_dist([lambda lambda],[1 2]);
yscale = sw_dist([1 2],lon([lon1 lon1]));

st = 5;
load invgray

% start figure
figure
colormap(invgray)

contourf(lon(lon1:lon2),lat(lat1:lat2),-Y(lat1:lat2,lon1:lon2),[0:250:2000],'linestyle','none')
caxis([0 3000])
set(gca,'plotboxaspectratio',[ratio 1 1])

pos_main = get(gca,'position');
h = cbarf(-Y(lat1:lat2,lon1:lon2),[0:250:2000]);
ylabel(h,'Depth (m)')
pos_cb = get(h,'position');
set(h,'position',[pos_cb(1) pos_main(2) pos_cb(3) pos_main(4)])

% plot mooring locations
hold on
plot(poslon,poslat,'k.','markersize',20)

for i = 1:nfile
    quiver(poslon(i),poslat(i),nanmean(unew(i,:),2)'*100/xscale,nanmean(vnew(i,:),2)'*100/yscale,'k','linewidth',1)
end

for di = 1:nfile
    [~,x,y] = varelip(uh(:,di),vh(:,di),0);
    plot(poslon(di)+100*x/xscale,poslat(di)+100*y/yscale,'k')
end

quiver(-22.8,67.6,10/xscale,0,'k')
text(-22.8,67.63,'10 cm/s')

ran = -21:21;
plot(poslon(3)+ran/xscale,poslat(3)+ran/sind(angB)/yscale,'k')

ran = 1;
% phase velocity
quiver(poslon(3),poslat(3),ran*cpkd*sind(mean(ang(2:4)))/xscale,ran*cpkd*cosd(mean(ang(2:4)))/yscale,'k','linewidth',2,'MaxHeadSize',10)
text(-23.68,67.6,'C_p','fontweight','bold')
quiver(poslon(3),poslat(3),ran*Cgkd*sind(Cgang)/xscale,ran*Cgkd*cosd(Cgang)/yscale,'k','linewidth',2,'linestyle','--')
text(-22.65,66.96,'C_g','fontweight','bold')

ox = [-.15 .09 .07 .07 -.33 -.33 -.33];
oy = [-.03 .01 .01 .01 .01 .01 .01];
for i = 1:length(poslon)
    text(poslon(i)+ox(i),poslat(i)+oy(i),['KGA-' num2str(i)])
end

% xlims = get(gca,'xlim');
% figposx = (poslon(3)-xlims(1))/diff(xlims);
% ylims = get(gca,'ylim');
% figposy = (poslat(3)-ylims(1))/diff(ylims);
% dx = ran*cpkd*sind(mean(ang(2:4)))/xscale;
% dy = ran*cpkd*cosd(mean(ang(2:4)))/yscale;
% annotation('arrow',figposx+[0 dx],figposy+[0 dy])%poslon(3)+[0 ran*Cgkd*sind(Cgang)/xscale],poslat(3)+[0 ran*Cgkd*cosd(Cgang)/yscale])

% plot(poslon(3)+ran/xscale,poslat(3)+ran/sind(mean(ang(2:4)))/yscale,'k')


axis([-25.6 -22.1 66.9 67.95]) 



print_fig('map_ellipse',savefold,1,.6)
