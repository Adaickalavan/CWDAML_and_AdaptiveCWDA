%Compute number of bit errors for 
%4-PSK,8-PSK,16-PSK,%8-circular,16-star,16-QAM,32-QAM,64-QAM,16-QAMGray,64-QAMGray 
%constellations  
%%
function bit_errors = count_error(format,est_sig,sig)

    err_index = find(est_sig ~= sig); %Find the symbol decisions which are different from the transmitted signals

    if strcmp(format,'PSK')
        error_bits = [1,1,2,1,2,2,3,1,2,2,3,2,3,3,4]; %Number of bits in error 
        gray = [0,1,3,2,6,7,5,4,12,13,15,14,10,11,9,8]; %Gray coded MPSK up to 16-PSK  
        current_error = error_bits(bitxor(gray(est_sig(err_index)+1), gray(sig(err_index)+1))); %Compute the number of bits in error
    elseif strcmp(format,'QAM') || ...
           strcmp(format,'QAMGray') || ...
           strcmp(format,'circular') || ...
           strcmp(format,'star') || ...
           strcmp(format,'triangular_NW') || ...
           strcmp(format,'triangular_mine')
        current_error = bitxor(est_sig(err_index), sig(err_index));  %Compute the number of bits in error
        current_error = sum(dec2bin(current_error)-48); %This line changes ASCII to numbers and 48 is ASCII value for '0'
    else
        error('Unavailable %s format',format)
    end

    bit_errors = sum(current_error);

end