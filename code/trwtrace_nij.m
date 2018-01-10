function [CG,WCHECK,LOC,K,BSLOPE,LAMBDA,H] = trwtrace_nij(loc,k,prd,parms)


ncorr=5;

% set up initial values
% prd = 4; %days
%loc = [-23.6 67.34];
 % k components cannot equal 0!
stepmax = 96;

minlon = -30;
maxlon = -18;
minlat = 65;
maxlat = 70;

% calculate splines
% load ~/data/WHOI/bathymetry/bathy
% bath = elev; lon = lons; lat = lats';
load ~/data/bathymetry/bathy_ibcao_ds


bath = bath(find(lat>minlat,1)-1:find(lat>maxlat,1),find(lon>minlon,1)-1:find(lon>maxlon,1));
lon = lon(find(lon>minlon,1)-1:find(lon>maxlon,1));
lat = lat(find(lat>minlat,1)-1:find(lat>maxlat,1));
[m,n] = size(bath);

for i = 1:n
    Y(:,i) = my_lowpass(bath(:,i),60);
end

for i = 1:m
    Y(i,:) = my_lowpass(Y(i,:),60);
end

bath = Y;
h = bath;
[h_x,h_y] = gradient(h,lon,lat);
[h_xx,h_xy] = gradient(h_x,lon,lat);
[h_yx,h_yy] = gradient(h_y,lon,lat);


% % assign parameters
% beta = parms(1);
% dt   = parms(2);
% f    = parms(3);
% H0   = parms(4);
% N0   = parms(5); 

% Set up parameters
beta = 10^-11; 
dt  = -.5; % in hours
f = 10^-4;
H0 = 1000;
N0 = 1e-3;

f = parms(1);
N0 = parms(2);
H0 = parms(3);

% nondimensionalize items and define useful conversion factors
RD=N0*H0/f;
fscale=RD*beta;
w=2*pi/(prd*24*3600*fscale);
dt=dt*3600*fscale;
fbet=f/fscale;
dgxy=[RD/86.263e+3 RD/111.0e+3]; xydg=dgxy;
dgxy=[1 dgxy dgxy.*dgxy dgxy(1)*dgxy(2)]./H0;
k=RD*k;

% initialize counters and accumulation arrays
clear K LOC CG WCHECK
step=0; K=k; LOC=loc; H=nan;
ncount=0;

% for now use uniform stratification
n=1; dndh=0;

while step<stepmax

% find topograpohy gradient
    hi = interp2(lon,lat,bath,loc(1),loc(2));
    h_xi = interp2(lon,lat,h_x,loc(1),loc(2));
    h_yi = interp2(lon,lat,h_y,loc(1),loc(2));
    h_xxi = interp2(lon,lat,h_xx,loc(1),loc(2));
    h_yyi = interp2(lon,lat,h_yy,loc(1),loc(2)); 
    h_xyi = interp2(lon,lat,h_xy,loc(1),loc(2));
    gradh = dgxy.*[hi h_xi h_yi h_xxi h_yyi h_xyi];
    
% save bottom slope
bslope=sqrt(gradh(2)^2+gradh(3)^2);

bg=fbet*gradh;
lam=n*sqrt(k(1)^2+k(2)^2+k(1)/w);
h=gradh(1); t=tanh(lam*h); g=lam*t; g1=t+lam*h*(1-t^2);
wcheck=(n^2*(bg(3)*k(1)-bg(2)*k(2))/(g*w))-1;
th=-(k(1)*n*n)/(2*w*w*lam); ph=g/(g1*w); cmn=1/(ph+th);
dwdh=((2*w*ph-lam)*dndh/n) - (lam^2*(1-t^2)/g1);
dwdhx=-(n^2*k(2)*ph)/g; dwdhy=-dwdhx*k(1)/k(2);


cg=cmn*[(ph*n^2*bg(3)/g)+(th*(w/k(1) + 2*w*w)),...
       (-ph*n^2*bg(2)/g)+(th*(2*k(2)*w*w/k(1)))];

dkdt=-cmn*[dwdh*gradh(2)+dwdhx*bg(4)+dwdhy*bg(6),...
      dwdh*gradh(3)+dwdhx*bg(6)+dwdhy*bg(5)];

step=step+1

if step~=1
ncount=ncount+1;
loc=LOC(step-1,:)+(xydg.*(2*dt*cg)); 
k=K(step-1,:)+2*dt*dkdt;

%trapezoidal correction
if ncount==ncorr
    % find topograpohy gradient
        hi = interp2(lon,lat,bath,loc(1),loc(2));
        h_xi = interp2(lon,lat,h_x,loc(1),loc(2));
        h_yi = interp2(lon,lat,h_y,loc(1),loc(2));
        h_xxi = interp2(lon,lat,h_xx,loc(1),loc(2));
        h_yyi = interp2(lon,lat,h_yy,loc(1),loc(2)); 
        h_xyi = interp2(lon,lat,h_xy,loc(1),loc(2));
        gradh = dgxy.*[hi h_xi h_yi h_xxi h_yyi h_xyi];
    bg=fbet*gradh;
    lam=n*sqrt(k(1)^2+k(2)^2+k(1)/w);
    h=gradh(1); t=tanh(lam*h); g=lam*t; g1=t+lam*h*(1-t^2);
    th=-(k(1)*n*n)/(2*w*w*lam); ph=g/(g1*w); cmn=1/(ph+th);
    dwdh=((2*w*ph-lam)*dndh/n) - (lam^2*(1-t^2)/g1);
    dwdhx=-(n^2*k(2)*ph)/g; dwdhy=-dwdhx*k(1)/k(2);

    cg_g=cmn*[(ph*n^2*bg(3)/g)+(th*(w/k(1) + 2*w*w)),...
           (-ph*n^2*bg(2)/g)+(th*(2*k(2)*w*w/k(1)))];

    dkdt_g=-cmn*[dwdh*gradh(2)+dwdhx*bg(4)+dwdhy*bg(6),...
              dwdh*gradh(3)+dwdhx*bg(6)+dwdhy*bg(5)];

    loc=LOC(step,:)+(xydg.*(dt*(cg+cg_g)/2));
    k=K(step,:)+dt*(dkdt+dkdt_g)/2;
    % find topograpohy gradient
        hi = interp2(lon,lat,bath,loc(1),loc(2));
        h_xi = interp2(lon,lat,h_x,loc(1),loc(2));
        h_yi = interp2(lon,lat,h_y,loc(1),loc(2));
        h_xxi = interp2(lon,lat,h_xx,loc(1),loc(2));
        h_yyi = interp2(lon,lat,h_yy,loc(1),loc(2)); 
        h_xyi = interp2(lon,lat,h_xy,loc(1),loc(2));
        gradh = dgxy.*[hi h_xi h_yi h_xxi h_yyi h_xyi];bg=fbet*gradh;
    lam=n*sqrt(k(1)^2+k(2)^2+k(1)/w);
    h=gradh(1); t=tanh(lam*h); g=lam*t;
    wcheck=(n^2*(bg(3)*k(1)-bg(2)*k(2))/(g*w))-1;
    ncount=0;
end

else
% euler time step
loc=LOC(step,:)+(xydg.*(dt*cg)); 
k=K(step,:)+dt*dkdt;

% trapezoidal correction
% find topograpohy gradient
    hi = interp2(lon,lat,bath,loc(1),loc(2));
    h_xi = interp2(lon,lat,h_x,loc(1),loc(2));
    h_yi = interp2(lon,lat,h_y,loc(1),loc(2));
    h_xxi = interp2(lon,lat,h_xx,loc(1),loc(2));
    h_yyi = interp2(lon,lat,h_yy,loc(1),loc(2)); 
    h_xyi = interp2(lon,lat,h_xy,loc(1),loc(2));
    gradh = dgxy.*[hi h_xi h_yi h_xxi h_yyi h_xyi];
bg=fbet*gradh;
lam=n*sqrt(k(1)^2+k(2)^2+k(1)/w);
h=gradh(1); t=tanh(lam*h); g=lam*t; g1=t+lam*h*(1-t^2);
th=-(k(1)*n*n)/(2*w*w*lam); ph=g/(g1*w); cmn=1/(ph+th);
dwdh=((2*w*ph-lam)*dndh/n) - (lam^2*(1-t^2)/g1);
dwdhx=-(n^2*k(2)*ph)/g;  dwdhy=-dwdhx*k(1)/k(2);

cg_g=cmn*[(ph*n^2*bg(3)/g)+(th*(w/k(1) + 2*w*w)),...
       (-ph*n^2*bg(2)/g)+(th*(2*k(2)*w*w/k(1)))];

dkdt_g=-cmn*[dwdh*gradh(2)+dwdhx*bg(4)+dwdhy*bg(6),...
          dwdh*gradh(3)+dwdhx*bg(6)+dwdhy*bg(5)];

loc=LOC(step,:)+(xydg.*(dt*(cg+cg_g)/2));
k=K(step,:)+dt*(dkdt+dkdt_g)/2;
% find topograpohy gradient
    hi = interp2(lon,lat,bath,loc(1),loc(2));
    h_xi = interp2(lon,lat,h_x,loc(1),loc(2));
    h_yi = interp2(lon,lat,h_y,loc(1),loc(2));
    h_xxi = interp2(lon,lat,h_xx,loc(1),loc(2));
    h_yyi = interp2(lon,lat,h_yy,loc(1),loc(2)); 
    h_xyi = interp2(lon,lat,h_xy,loc(1),loc(2));
    gradh = dgxy.*[hi h_xi h_yi h_xxi h_yyi h_xyi];
bg=fbet*gradh;
lam=n*sqrt(k(1)^2+k(2)^2+k(1)/w);
h=gradh(1); t=tanh(lam*h); g=lam*t;
wcheck=(n^2*(bg(3)*k(1)-bg(2)*k(2))/(g*w))-1;
end

if step~=1
CG=[CG;sqrt(cg(1)^2+cg(2)^2)];
WCHECK=[WCHECK;wcheck];
BSLOPE=[BSLOPE;bslope];
LAMBDA=[LAMBDA;lam];

else
CG=sqrt(cg(1)^2+cg(2)^2);
WCHECK=wcheck;
BSLOPE=bslope;
LAMBDA=lam;
end

LOC=[LOC;loc];
K=[K;k];
H = [H;h];

% Put in a cutoff so that the trace cannot continue into an area where
% it is not deep enough for the assumptions to remain true.  Note that 
% the depth is non-dimensionalized here by a factor of H.
% if h > 1
%    step=stepmax;
% end

% Put in cuttoffs so that trace cannot go out of spline area.
if loc(1) > maxlon
  step=stepmax;
end
if loc(1) < minlon
  step=stepmax;
end
if loc(2) > maxlat
  step=stepmax;
end
if loc(2) < minlat
  step=stepmax;
end

end % while
CG=(RD*fscale*3.6*24)*CG;

% figure
% contourf(lon,lat,bath,[-4000:1000:-1000 -500:100:500],'linestyle','none')
% hold on
% contour(lon,lat,bath,[0 0],'k')
% plot(LOC(:,1),LOC(:,2),'k')
