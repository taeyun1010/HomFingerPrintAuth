function [output] = reducedoublebits(input, intbit, fracbit)
%REDUCEDOUBLEBITS 

bits = f_d2b(input);

strn=strtrim(num2str(bits));
     
%------------------------------------------------
if isempty(find(strn=='.'))  
    if length(bits) > intbit
        reduceddouble=bits(length(bits)-intbit+1:end);
    else
        reduceddouble = bits;
    end      
else
%------------------------------------------------
    k = find(strn =='.');

    %Retrieving INTEGER and FRACTIONAL PARTS as strings
    i_part=strn(1:k-1);
    f_part=strn(k:end);

    %reduce the number of bits to intbit and fracbit
    if length(i_part) > intbit
        i_part=i_part(length(i_part)-intbit+1:end);
    end
    if length(f_part) > fracbit
        f_part=f_part(1:fracbit+1);
    end
    % %Converting the strings back to numbers
    % ni_part=str2num(i_part);
    % nf_part=str2num(f_part);

    reduceddouble = strcat(i_part, f_part);
end
    output = f_b2d(reduceddouble);


end

