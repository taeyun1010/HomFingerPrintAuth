function [output] = addtion(encryptedBits1, encryptedBits2)
%ADDTION Summary of this function goes here
%   Detailed explanation goes here

% numberof rows
numBits1 = size(encryptedBits1,1);
numBits2 = size(encryptedBits2,1);

%initialize carry bit to zero and encrypt
[carrya, carryb] = encMat(0);
encryptedzeroa = carrya;
encryptedzerob = carryb;

if numBits1 < numBits2
    
    for i = 1:numBits1
        a1 = encryptedBits1{i,1};
        %     disp(a{1});
        %     disp(a{2});
        %     disp(a{3});
        %temp = a(1);
        %     for j = 1:500
        %         disp(a(j));
        %     end
        b1 = encryptedBits1{i,2};
        a2 = encryptedBits2{i,1};
        b2 = encryptedBits2{i,2};
        
        [suma,sumb, carrya, carryb] = addBit(a1,b1,a2,b2,carrya,carryb);
        
        output{i,1} = suma;
        output{i,2} = sumb;
        
    end
    
    
    
    % add to remaining carry
    for i = (numBits1 + 1):numBits2
        if i == (numBits1 + 1)
            % add 0 bit as a first argument for the first loop
            a1 = encryptedzeroa;
            b1 = encryptedzerob;
        else
             a1 = encryptedBits1{i,1};
             b1 = encryptedBits1{i,2};
        end
        a2 = encryptedBits2{i,1};
        b2 = encryptedBits2{i,2};
        [suma, sumb, carrya, carryb] = addBit(a1,b1,a2,b2,carrya,carryb);
        output{i,1} = suma;
        output{i,2} = sumb;
        
        % if last loop, must pass carry (if any) to the upper bit
        if i==numBits2
            output{(i+1),1} = carrya;
            output{(i+1),2} = carryb;
        end
        
    end
    
elseif numBits1 == numBits2
      for i = 1:numBits1
        a1 = encryptedBits1{i,1};
        %     disp(a{1});
        %     disp(a{2});
        %     disp(a{3});
        %temp = a(1);
        %     for j = 1:500
        %         disp(a(j));
        %     end
        b1 = encryptedBits1{i,2};
        a2 = encryptedBits2{i,1};
        b2 = encryptedBits2{i,2};
        
        [suma,sumb, carrya, carryb] = addBit(a1,b1,a2,b2,carrya,carryb);
        
        output{i,1} = suma;
        output{i,2} = sumb;
        
      end
      
      % if last loop, must pass carry (if any) to the upper bit
      if i==numBits1
          output{(i+1),1} = carrya;
          output{(i+1),2} = carryb;
      end
    
      
      % input1 has more bits
else
    
    for i = 1:numBits2
        a1 = encryptedBits1{i,1};
        %     disp(a{1});
        %     disp(a{2});
        %     disp(a{3});
        %temp = a(1);
        %     for j = 1:500
        %         disp(a(j));
        %     end
        b1 = encryptedBits1{i,2};
        a2 = encryptedBits2{i,1};
        b2 = encryptedBits2{i,2};
        
        [suma,sumb, carrya, carryb] = addBit(a1,b1,a2,b2,carrya,carryb);
        
        output{i,1} = suma;
        output{i,2} = sumb;
        
    end
    
    
    
    % add to remaining carry
    for i = (numBits2 + 1):numBits1
        if i == (numBits2 + 1)
            % add 0 bit as a first argument for the first loop
            a2 = encryptedzeroa;
            b2 = encryptedzerob;
        else
             a2 = encryptedBits1{i,1};
             b2 = encryptedBits1{i,2};
        end
        a1 = encryptedBits1{i,1};
        b1 = encryptedBits1{i,2};
        [suma, sumb, carrya, carryb] = addBit(a1,b1,a2,b2,carrya,carryb);
        output{i,1} = suma;
        output{i,2} = sumb;
        
        % if last loop, must pass carry (if any) to the upper bit
        if i==numBits1
            output{(i+1),1} = carrya;
            output{(i+1),2} = carryb;
        end
        
    end
    
end

end

