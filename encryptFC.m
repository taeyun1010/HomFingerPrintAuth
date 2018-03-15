function partialoutput = encryptFC(fingercode)
%ENCRYPTFC Summary of this function goes here
%   Detailed explanation goes here

% n = number of disks
[m,n] = size(fingercode);

for i = 1:n
    partialfc = fingercode{1,i};
    
    % s = number of sectors
    [s,t] = size(partialfc);
    for j = 1:s
        encryptedBits = homEncrypt(round(partialfc(j)));
        partialoutput{1,j} = encryptedBits;
    end
    
end

end

