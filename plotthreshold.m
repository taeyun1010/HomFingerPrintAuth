maxthreshold = max(max(distances), max(distances2));

thresholds = [0:1:maxthreshold];
fns = zeros(1+maxthreshold);
tps = zeros(1+maxthreshold);
fps = zeros(1+maxthreshold);
tns = zeros(1+maxthreshold);
for i=1:(1+maxthreshold)
     [fn,tp, fp, tn] = falsenegative(distances,distances2,thresholds(i));
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