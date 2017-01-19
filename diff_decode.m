%Differential decoding of 
%8-circular,16-circular,16-star,16-QAMGray,16-QAM,32-QAM,64-QAM,64-QAMGray
%constellations
%%
function [est_sig, prev_rx_sector] = diff_decode(index,prev_rx_sector,S,sector_rotation,const_points)

    est_sig = zeros(length(index),1);
    for i = 1:length(index)
        cur_rx_sector = sector_rotation(index(i));
        decoded_rx_data = const_points(index(i))*cur_rx_sector/prev_rx_sector;
        prev_rx_sector = cur_rx_sector; %Reset previous received constellation point
        [unwanted,est_sig(i,1)] = min(abs(decoded_rx_data - S));
    end
    %To ensure that est_sig is in the range of [0, M-1]
    est_sig = est_sig - 1;

end