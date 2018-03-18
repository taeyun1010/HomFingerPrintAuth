maxthreshold = max(max(distances4), max(distances3));
maxthreshold = int32(maxthreshold);
thresholds = [0:1:maxthreshold];
fns = zeros(1+maxthreshold);
tps = zeros(1+maxthreshold);
fps = zeros(1+maxthreshold);
tns = zeros(1+maxthreshold);
for i=1:(1+maxthreshold)
     [fn,tp, fp, tn] = falsenegative(distances4,distances3,thresholds(i));
     fns(i) = fn;
     tps(i) = tp;
     fps(i) = fp;
     tns(i) = tn;
end

figure
plot(thresholds, fns, 'c');
hold on;
 plot(thresholds, tps, 'r');
 plot(thresholds, fps, 'k');
 plot(thresholds, tns, 'm');

xlabel('Threshold');
ylabel('fn,tp,fp,tn');
legend('False Negative','True Positive','False Positive','True Negative');