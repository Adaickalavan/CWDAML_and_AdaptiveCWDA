%Differential encoding of 
%8-circular,16-circular,16-star,16-QAMGray,16-QAM,32-QAM,64-QAM,64-QAMGray
%constellations
%%
function [tx_data, prev_tx_sector] = diff_encode(sig,prev_tx_sector,sector_rotation,const_points) 

    tx_data = zeros(length(sig),1);
    for i = 1:length(sig)
        sig_sector = sector_rotation(sig(i)+1);
        cur_tx_sector = sig_sector*prev_tx_sector;
        tx_data(i) = const_points(sig(i)+1)*cur_tx_sector;
        prev_tx_sector = cur_tx_sector; 
    end

end