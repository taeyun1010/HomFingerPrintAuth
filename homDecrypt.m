function [plaintext] = homDecrypt(encryptedBits)
%HOMDECRYPT Summary of this function goes here
%   Detailed explanation goes here

% numberof rows
numBits = size(encryptedBits,1);
plaintext = 0;


for i = 1:numBits
    a = encryptedBits{i,1};
%     disp(a{1});
%     disp(a{2});
%     disp(a{3});
    %temp = a(1);
    %     for j = 1:500
%         disp(a(j));
%     end
    b = encryptedBits{i,2};
    bit = decMat(a,b);
    plaintext = plaintext + bit * 2^(i-1);
    
end


end

