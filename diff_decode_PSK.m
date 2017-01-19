%Differentially decode M-PSK
%%
function [est_sig, rx_pre_sig] = diff_decode_PSK(M,rx_cur_sig,rx_pre_sig)

    decoded_sig = rx_cur_sig - [rx_pre_sig; rx_cur_sig(1:end-1)]; %Differentially decode the signal
    rx_pre_sig = rx_cur_sig(end); %Reset previous received constellation point          
    est_sig = mod(decoded_sig, M); %Differential decoding        

end