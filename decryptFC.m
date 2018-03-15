function fingercode = decryptFC(encryptedFC)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

% n = number of disks
[m,n] = size(encryptedFC);

for i = 1:n
    partialEfc = encryptedFC{1,i};
    
    % s = number of sectors
    [s,t] = size(partialEfc);
    for j = 1:s
        encryptedBits = partialEfc{j,1};
        partialoutput{j,1} = homDecrypt(encryptedBits);
    end
    fingercode{1,i} = partialoutput;
end

end

