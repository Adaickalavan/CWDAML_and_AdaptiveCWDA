%Karpaga Vinayagar - Pillaiyarpatti
%Adaickalavan Meiyappan - NUS

%%
function [ave_BER] = NL_Adaptive_CWDA_MPSK(rp)

disp('--------------------');
disp('Start of program run');
format shortG %Number display format
time1 = clock; %Retrieve starting time to compute code run time later

%Set the random number stream
reset(rp.stream);

%%
M = rp.M; % QAM/QAMGray/circular/star
bit_no = log2(M); %Number of bits in each symbol
S = constellation(rp); %Generate the signal constellation 

training_length = rp.training_length; %Number of training symbols
R = rp.bit_rate/bit_no; %Symbol rate, Units: symbols/s
laser_linewidth = rp.laser_linewidth; %Combined laser linewdith, Units: Hz
rp.sigma = sqrt(2*pi*laser_linewidth/R); %Standard deviation of laser phase noise
frequency_offset = rp.frequency_offset; %Constant frequency offset, Units: Hz
rp.omega = 2*pi*frequency_offset/R; %Constant phase offset due to frequency offset, Units: rad
% rp.omega_FM = 2*pi*rp.freq_mod/R; %Frequency modulation due to frequency fluctuation, Units: rad
total_run = rp.total_run; %Number of times to repeat same simulation. Mainly to obtain ensemble averages.

rp = NL_par(rp); % Compute nonlinear parameter list

fprintf('LLW = %4.3e\n',laser_linewidth);
fprintf('FO = %4.3e\n',frequency_offset);
fprintf('Bo = %4.3e\n',rp.Bo);
% SNR_bit_dB = 10*log10((1e-3*10.^(rp.power/10))./(rp.NA*rp.P_ASE)) - 10*log10(log2(rp.M));
% fprintf('SNR_bit = %4.2f dB\n\n',SNR_bit_dB);

%%

BER = zeros(length(rp.power),total_run); %Preset BER matrix for faster computation

%Run the code several times to obtain ensemble average values for BEP
for run = 1:total_run
    
    power_index = 0; %Initialize power index pointer

    %Repeat until Bit Error Probability is below 1e-5                      
    while power_index < length(rp.power)  
        
        %Loop paramaters
        power_index = power_index + 1; %Increment index pointer
        power = rp.power(power_index); %Increase the power by step size
        A = sqrt(1e-3*10^(power/10)); %Signal amplitude, Power_linear = |A|^2
        
        %Initialization
        tx_pre_sig = 0; %Initial constellation point sent
        rx_pre_sig = tx_pre_sig; %Previously sent constellation point
        phase_noise = 0; %Initial cumulative phase noise
        input_cur = zeros(2,1); %Initialize x(k)
        r_vector = zeros(2,1);
        R_matrix = 0.01*eye(2);
        w_cur = [0; 1]; %Initialize w(k)
        V = 1;
        No_sent_symbols = 0; %Initialize total number of sent bits to zero
        total_error_bit = 0; %Initialize total bits received in error to zero
        
        %Repeat signal transmission and reception until error bits >= bit_no*100
        while No_sent_symbols < rp.No_sent_symbols && total_error_bit < rp.total_error_bit
            
            No_sent_symbols = No_sent_symbols + 1; %Increment the total sent bits
            %Generate a signal point to transmit
            sig = randi([0,M-1]); %Generate a constellation point to send
            %Differential encoding
            [tx_cur_sig, tx_pre_sig] = diff_encode_PSK(M,sig,tx_pre_sig); %Diferential encoding

            %Received message
            m = A*S(tx_cur_sig + 1); %Actual sent message    
%             rp.noss = No_sent_symbols;
            [rx, phase_noise] = AWGN_NLPN_LPN_FO_channel(m,phase_noise,rp); %AWGN, NLPN, LPN, FO impairment
                  
            %Training period
            if No_sent_symbols <= training_length
                rx_pre_sig = tx_pre_sig; %Replace estimated received constellation point
                m_cap = m; %Replace estimated message with actual message 
            else
                %Decision aided maximum likelihood coherent symbol by symbol detection
                dist = abs(rx*conj(V)/A - S(1:M)); %Compute distance between detected point and constellations point
                [~,index] = min(dist);
                rx_cur_sig = index - 1; %Current received constellation point
                m_cap = A*S(index); %Decision on received message              
            end
            %Ideal decision feedback
            if rp.idl_dcs_fdb == 1
                m_cap = m; 
            end
            
            %Form filter input vector 
            input_prev1 = input_cur; %Form x(k-1)
            input_cur = [V; rx/m_cap]; %Form x(k)

            %Tap weight adaptation
            if No_sent_symbols > 1 && No_sent_symbols <= rp.freq_est_length
                r_vector = r_vector + (rx/m_cap)*conj(input_prev1); 
                R_matrix = R_matrix + conj(input_prev1)*(input_prev1.');
                w_cur = R_matrix\r_vector; %w is a column vector (K rows) of LMMSE filter coefficients
            end
            
            %Form reference phasor
            V = (w_cur.')*input_cur; %Form V(k+1)
            
            %Decode and compute bit error
            if No_sent_symbols > training_length
                %Differential decoding
                [est_sig, rx_pre_sig] = diff_decode_PSK(M,rx_cur_sig,rx_pre_sig); %Differential decoding
                %Error counting
                bit_errors = count_error(rp.format,est_sig,sig);
                total_error_bit = total_error_bit + bit_errors;           
            end

        end      
        
        fprintf('%s\n',datestr(now));
        BEP = total_error_bit/(bit_no*(No_sent_symbols-training_length)); %Compute bit error probability
        fprintf('Power = %3.1f dBm\n',power);
        fprintf('Sent symbols = %u\n',No_sent_symbols);
        fprintf('Bit errors = %u\n',total_error_bit);
        fprintf('BEP = %7.5e\n\n',BEP);
        BER(power_index,run) = BEP; %Store the BEP  

    end
    
end

%%
%Compute the ensemble average of BER vs SNR per bit
ave_BER = mean(BER,2); %Averaging along each row  

%Compute SNR value at specified BER value
y0 = ave_BER;
x0 = rp.power.';
power_value = interpolate_NL(rp,x0,y0);
fprintf('Fall-1 @ %6.4f dBm\n', power_value); %Crossing at falling edge of rp.BEP(1)
% fprintf('Rise-1 @ %6.4f dBm\n', power_value(2,1)); %Crossing at rising edge of rp.BEP(1)
% fprintf('Fall-2 @ %6.4f dBm\n', power_value(1,2)); %Crossing at falling edge of rp.BEP(2)
% fprintf('Rise-2 @ %6.4f dBm\n', power_value(2,2)); %Crossing at rising edge of rp.BEP(2)

%%
%Save all variables from current workspace for later access
full_name = mfilename('fullpath'); %Obtain the full name of this script/function 
[~, filename, ~] = fileparts(full_name); %Parse out the filename alone
pathname = './Results/'; %Set the location path to save
save_loc = [pathname,filename,'_',num2str(M,'%u'),'-',rp.format,'_BR_',num2str(R,'%10.3G'),'_LLW_',num2str(laser_linewidth,'%10.3G'),'_FO_',num2str(frequency_offset,'%10.3G'),'_NA_',num2str(rp.NA,'%u'),'.mat'];
save(save_loc);

%%
%Get the screensize to specify figure size and location
scrsz = get(0,'ScreenSize'); 

% Specify position of figure on screen. rect = [left, bottom, width, height]
figure('OuterPosition',[1 scrsz(4)/2 scrsz(3)/4 scrsz(4)/2],'Name','BER vs SNR_b') 
%Plot the BER vs SNR graph
semilogy(rp.power,ave_BER,'-b.','linewidth',1,'markerfacecolor','r')
title('BER vs Launch power'),
xlabel('Launch power (dBm)'),ylabel('BER'),grid;    

% % Specify position of figure on screen. rect = [left, bottom, width, height]
% figure('OuterPosition',[1 scrsz(4)/2 scrsz(3)/4 scrsz(4)/2],'Name','BER vs SNR_b') 
% %Plot the BER vs SNR graph
% semilogy(SNR_bit_dB,ave_BER,'-b.','linewidth',1,'markerfacecolor','r')
% title('BER vs SNR per bit(dB)'),
% xlabel('SNR per bit (dB)'),ylabel('BER'),grid;  

%%
%Find the elapsed time for the code to run
time2 = clock;
elapsed_minutes = etime(time2, time1)/60;
fprintf('Elapsed Time = %6.2f minutes\n', elapsed_minutes);
disp('End of program run');
disp('------------------');

end
