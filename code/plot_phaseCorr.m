figure
col = {'k','r','g','b','c','y'};
for i = 1:6
    [varcorr,lags] = myxcorr(wsr(:,i),wsr(:,i+1),0,0,1/8);
    varcorr = varcorr(end:-1:1);
    lags = lags(end:-1:1);
    midi = (length(varcorr)+1)/2;
    maxlag(i) = lags(midi - 1 + find(varcorr(midi:end)==max(varcorr(midi:end))));
    plot(lags,varcorr,col{i})
    hold on
    
    axis([-15 15 -1 1])
end
if i == 6
    legend('1-2','2-3','3-4','4-5','5-6','6-7')
else
    legend('1-2','2-3','3-4','4-5')
end 
plot([lags(1) lags(end)],[0 0],'k')
plot([0 0],[-1 1],'k')

