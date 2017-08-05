%Karpaga Vinayagar - Pillaiyarpatti
%Adaickalavan Meiyappan - NUS

%%
function [SNR_value] = Adaptive_CWDA_MPSK(rp)

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
sigma = sqrt(2*pi*laser_linewidth/R); %Standard deviation of laser phase noise
frequency_offset = rp.frequency_offset; %Constant frequency offset, Units: Hz
omega = 2*pi*frequency_offset/R; %Constant phase offset due to frequency offset, Units: rad
SNR_begin = rp.SNR_begin; %Signal to noise ratio per bit, Units: dB
SNR_max = rp.SNR_max; %Maximum SNR limit before declaring a phase tracking failure
step_size = rp.step_size; %Incremental value of SNR, Units: dB
total_run = rp.total_run; %Number of times to repeat same simulation. Mainly to obtain ensemble averages.

fprintf('LLW = %4.3e\n',laser_linewidth);
fprintf('FO = %4.3e\n\n',frequency_offset);

%%

highest_SNR = 0; %Highest achieved SNR over all runs
BER = zeros((1/step_size)*(SNR_max-SNR_begin)+1,total_run); %Preset BER matrix for faster computation

%Run the code several times to obtain ensemble average values for BEP
for run = 1:total_run
    
    BEP = 1; %Initialize bit error probability to 1
    SNR = SNR_begin - step_size; %Initialize SNR to SNR_begin - 0.5dB
    SNR_index = 0; %Initialize SNR index pointer
    
    %Repeat until Bit Error Probability is below 1e-5                      
    while BEP >= rp.BEP   
        
        %Loop paramaters
        SNR = SNR + step_size; %Increase the SNR by 0.5dB 
        SNR_index = SNR_index + 1; %Increment index pointer
    %   Es/(bit_no*No) = SNR_per_bit
    %   Es = bit_no*(10^(SNR_per_bit_in_dB/10)), where No = 1
        A = sqrt(bit_no*(10^(SNR/10))); %Signal amplitude, A = sqrt(Es)
        
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
            phase_noise = sigma*randn + phase_noise + omega;  %Current cumulative phase noise, Units: radians
            rx = m*exp(1j*phase_noise) + (randn + 1j*randn)/sqrt(2); %Received signal with laser phase noise and AWGN
            %Var['randn'] = 1  and Var['randn'/sqrt(2)] = 1/2 = No/2 as No = 1            
                                          
            %Training period
            if No_sent_symbols <= training_length
                rx_pre_sig = tx_pre_sig; %Replace estimated received constellation point
                m_cap = m; %Replace estimated message with actual message 
            else
                %Decision aided maximum likelihood coherent symbol by symbol detection
                dist = abs(rx*conj(V)/A - S(1:M)); %Compute distance between detected point and constellations point
                [unwanted,index] = min(dist);
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
        fprintf('SNR per bit = %3.1f dB\n',SNR);
        fprintf('Sent symbols = %u\n',No_sent_symbols);
        fprintf('Bit errors = %u\n',total_error_bit);
        fprintf('BEP = %7.5e\n\n',BEP);
        BER(SNR_index,run) = BEP; %Store the BEP  
        
        %Exit loop if the tracking has failed. Prevent infinite loop.
        if SNR >= SNR_max; 
            break
        end
        
    end
    
    %Update the highest SNR achieved over all runs
    if highest_SNR < SNR
        highest_SNR = SNR;
    end
    
end

%%
%Compute the ensemble average of BER vs SNR per bit
ave_BER = mean(BER,2); %Averaging along each row  

%Compute SNR value at specified BER value
y0 = ave_BER(1:(1/step_size)*(highest_SNR-SNR_begin)+1);
x0 = (SNR_begin:step_size:highest_SNR).';
% SNR_value = interpolate(rp,x0,y0);
SNR_value = 0;
fprintf('SNR at BER = %5.2e is %6.4f dB\n', rp.read_BEP, SNR_value);

%%
%Save all variables from current workspace for later access
full_name = mfilename('fullpath'); %Obtain the full name of this script/function 
[unwanted, filename, unwanted] = fileparts(full_name); %Parse out the filename alone
pathname = './Simulation Results/'; %Set the location path to save
save_loc = [pathname,filename,'_',num2str(M,'%u'),'-PSK_BR_',num2str(R,'%10.3G'),'_LLW_',num2str(laser_linewidth,'%10.3G'),'_FO_',num2str(frequency_offset,'%10.3G'),'_idlFdb_',num2str(rp.idl_dcs_fdb,'%u'),'.mat'];
save(save_loc);

%%
%Get the screensize to specify figure size and location
scrsz = get(0,'ScreenSize'); 

% Specify position of figure on screen. rect = [left, bottom, width, height]
figure('OuterPosition',[1 scrsz(4)/2 scrsz(3)/4 scrsz(4)/2],'Name','BER vs SNR_b') 
%Plot the BER vs SNR graph
semilogy(SNR_begin:step_size:highest_SNR,ave_BER(1:(1/step_size)*(highest_SNR-SNR_begin)+1),'-r.','linewidth',1,'markerfacecolor','r')
title('BER against SNR'),
xlabel('SNR_b (E_b/N_0) dB'),ylabel('BER'),grid;    

%%
%Find the elapsed time for the code to run
time2 = clock;
elapsed_minutes = etime(time2, time1)/60;
fprintf('Elapsed Time = %6.2f minutes\n', elapsed_minutes);
disp('End of program run');
disp('------------------');

end
