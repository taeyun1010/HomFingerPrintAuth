function [suma,sumb,carrya,carryb] = addBit(input1a,input1b, input2a, input2b, inputcarrya, inputcarryb)
%ADDBIT Summary of this function goes here
%   Detailed explanation goes here
    
% input1 XOR input2
[1xor2a, 1xor2b] = xorMat(input1a, input1b, input2a, input2b);

% sum, without carry
[suma, sumb] = xorMat(1xor2a, 1xor2b, inputcarrya, inputcarryb);

% input1 AND input2
[1and2a, 1and2b] = andMat(input1a, input1b, input2a, input2b);

% (input1 xor input2) and inputcarry
[intermediatea, intermediateb] = andMat(1xor2a, 1xor2b, inputcarrya, inputcarryb);

% carry
[carrya, carryb] = orMat(1and2a, 1and2b, intermediatea, intermediateb);

end

