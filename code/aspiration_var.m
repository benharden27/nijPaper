clear
load ~/data/KGA/products/all_gridded.mat

D = repmat(dvec,[length(tvec) 1 length(xvec)]);
X = permute(repmat(xvec,[length(tvec) 1 length(dvec)]),[1 3 2]);
aspi = D>650 & X<88;

wsasp = ws;
wsasp(~aspi) = nan;
wcasp = wc;
wcasp(~aspi) = nan;

egct = ws(:,1:13,7:11);
egct = nanmean(egct(:,:),2);
egct_c = wc(:,1:13,7:11);
egct_c = nanmean(egct_c(:,:),2);


egcb = ws(:,13:end,7:11);
egcb = nanmean(egcb(:,:),2);
egcb_c = wc(:,13:end,7:11);
egcb_c = nanmean(egcb_c(:,:),2);






dx = diff(xvec([1 2]));
dd = diff(dvec([1 2]));

nlocs = length(find(~isnan(nanmean(wsasp(:,:)))));

wsasp = nanmean(wsasp(:,:),2);
wcasp = nanmean(wcasp(:,:),2);



figure, 
subplot(2,1,1), hold on
plot(egct)
plot(egcb,'r')

subplot(2,1,2), hold on
plot(egct_c)
plot(egcb_c,'r')



figure, hold on
% plot(wsasp)
plot(wcasp,'r')