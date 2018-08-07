maxthreshold = max(max(distances10), max(distances9));
maxthreshold = int32(maxthreshold);
thresholds = [0:100000:maxthreshold];
temp = size(thresholds);
numberofthresholds = temp(2);
fns = zeros(numberofthresholds);
tps = zeros(numberofthresholds);
fps = zeros(numberofthresholds);
tns = zeros(numberofthresholds);



for i=1:(numberofthresholds)
     [fn,tp, fp, tn] = falsenegative(distances10,distances9,thresholds(i));
     fns(i) = fn;
     tps(i) = tp;
     fps(i) = fp;
     tns(i) = tn;
end

figure
plot(thresholds, fns, 'c');
hold on
 plot(thresholds, tps, 'r');
 plot(thresholds, fps, 'k');
 plot(thresholds, tns, 'g');

xlabel('Threshold');
ylabel('fn,tp,fp,tn');
legend('False Negative','True Positive','False Positive','True Negative' );
