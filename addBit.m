function [suma,sumb,carrya,carryb] = addBit(input1a,input1b, input2a, input2b, inputcarrya, inputcarryb)
%ADDBIT Summary of this function goes here
%   Detailed explanation goes here
    
%tic

% input1 XOR input2
[onexor2a, onexor2b] = xorMat(input1a, input1b, input2a, input2b);

% sum, without carry
[suma, sumb] = xorMat(onexor2a, onexor2b, inputcarrya, inputcarryb);

% input1 AND input2
[oneand2a, oneand2b] = andMat(input1a, input1b, input2a, input2b);

% (input1 xor input2) and inputcarry
[intermediatea, intermediateb] = andMat(onexor2a, onexor2b, inputcarrya, inputcarryb);

% carry
[carrya, carryb] = orMat(oneand2a, oneand2b, intermediatea, intermediateb);

%toc

end

