function [numFalseNegative,numTruePositive, numFalsePositive, numTrueNegative] = falsenegative(distances,distances2,threshold)
%FALSENEGATIVE distances must contain distances of fps of same person
%   Detailed explanation goes here

numFalseNegative = 0;
numTruePositive = 0;
numFalsePositive = 0;
numTrueNegative = 0;

% TODO : CAUTION: hard coded length
for i=1:2000
    distance = distances(i);
    if distance > threshold
        numFalseNegative = numFalseNegative + 1;
    else
        numTruePositive = numTruePositive + 1;
        
    end
end

% TODO: CAUTION
for i=1:1999
    distance = distances2(i);
    if distance > threshold
        numTrueNegative = numTrueNegative + 1;
    else
        numFalsePositive = numFalsePositive + 1;
    end
end

end

