%Differentially encode M-PSK
%%
function [tx_cur_sig, tx_pre_sig] = diff_encode_PSK(M,sig,tx_pre_sig)

    tx_cur_sig = mod(tx_pre_sig + cumsum(sig), M); %Diferential encoding
    tx_pre_sig = tx_cur_sig(end); %Reset previous sent constellation point      
              
end