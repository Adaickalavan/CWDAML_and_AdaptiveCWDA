%Karpaga Vinayagar - Pillaiyarpatti
%Adaickalavan Meiyappan - NUS

%%
function [ave_BER] = NL_CW_recursiveW_DAML_MPSK(rp)

disp('--------------------');
disp('Start of program run');
format shortG %Number display format
time1 = clock; %Retrieve starting time to compute code run time later

%Set the random number stream
reset(rp.stream);

%%
M = rp.M; % 16-QAM
bit_no = log2(M); %Number of bits in each symbol
S = constellation(rp);    

N = rp.filter_length; %Filter length
training_length = rp.training_length; %Number of training symbols
R = rp.bit_rate/bit_no; %Symbol rate, Units: symbols/s
laser_linewidth = rp.laser_linewidth; %Combined laser linewdith, Units: Hz
rp.sigma = sqrt(2*pi*laser_linewidth/R); %Standard deviation of laser phase noise
frequency_offset = rp.frequency_offset; %Constant frequency offset, Units: Hz
rp.omega = 2*pi*frequency_offset/R; %Constant phase offset due to frequency offset, Units: rad
total_run = rp.total_run; %Number of times to repeat same simulation, to obtain ensemble averages.

rp = NL_par(rp); % Compute nonlinear parameter list

fprintf('LLW = %4.3e\n',laser_linewidth);
fprintf('FO = %4.3e\n',frequency_offset);
fprintf('Bo = %4.3e\n',rp.Bo);
fprintf('FL = %u\n',N);
% SNR_bit_dB = 10*log10((1e-3*10.^(rp.power/10))./(rp.NA*rp.P_ASE)) - 10*log10(log2(rp.M));
% fprintf('SNR_bit = %4.2f dB\n\n',SNR_bit_dB);

%%
if training_length < 2*N
    error('WARNING: Training length is shorter than filter length N')
end

BER = zeros(length(rp.power),total_run); %Preset BER matrix for faster computation
% store_max = 5e3;
% store_rx = zeros(store_max,1);
% store_sig = zeros(store_max,1);

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
        phase_noise = 0; %Initial cumulative phase noise
        tx_pre_sig = 0; %Initial constellation point sent
        rx_pre_sig = tx_pre_sig; %Previously sent constellation point
        input_cur = zeros(N,1); %Initialize x(k)
        rx_m_cap = zeros(N,1);
%         r_vector = zeros(N,1);
%         R_matrix = 0.01*eye(N);
        R_inv = (1/0.001)*eye(N); %Inverse autocorrelation matrix
        w_cur = [1; zeros(N-1,1)]; %Initialize w(k)
        No_sent_symbols = 0; %Initialize total number of sent bits to zero
        total_error_bit = 0; %Initialize total bits received in error to zero
               
        %Repeat signal transmission and reception until error bits >= bit_no*100
        while No_sent_symbols < rp.No_sent_symbols && total_error_bit < rp.total_error_bit
            
            No_sent_symbols = No_sent_symbols + 1; %Increment the total sent bits
            %Generate a signal point to transmit
            sig = randi([0,M-1]); %Generate a constellation point to send
            %Differential encoding of data for 16-QAM
            [tx_cur_sig, tx_pre_sig] = diff_encode_PSK(M,sig,tx_pre_sig); %Diferential encoding

            %Received message
            m = A*S(tx_cur_sig + 1); %Actual sent message
            [rx, phase_noise] = AWGN_NLPN_LPN_FO_channel(m,phase_noise,rp); %AWGN, NLPN, LPN, FO impairment

%             %Store received signals for graph plotting
%             if run == 1 && power == power_max && No_sent_symbols <= store_max
%                 store_rx(No_sent_symbols,1) = rx;
%                 [~,I_temp] = min(abs(tx_data - S));
%                 store_sig(No_sent_symbols,1) = I_temp;
%             end 
            
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
                    
            %Form filter input vector 
            input_prev1 = input_cur; %Form x(k-1)
            input_cur = [rx*conj(m_cap); input_cur(1:end-1)]; %Form x(k)
            rx_m_cap = [m_cap; rx_m_cap(1:end-1)];

            %Tap weight adaptation
            if No_sent_symbols > 1 && No_sent_symbols <= rp.freq_est_length
%                 r_vector = r_vector + c*(rx/m_cap)*conj(input_prev1); 
%                 R_matrix = R_matrix + (c^2)*conj(input_prev1)*(input_prev1.');
%                 w_cur = R_matrix\r_vector; %w is a column vector (K rows) of LMMSE filter coefficients
            
                %Recursive weight vector computation
                intermediate_vector = c*R_inv*conj(input_prev1);
                gain_vector = intermediate_vector./(1 + c*(input_prev1.')*intermediate_vector);
                error_estimate = rx/m_cap - c*(w_prev1.')*input_prev1;
                w_cur = w_prev1 + gain_vector.*error_estimate;
                R_inv = R_inv - gain_vector*(intermediate_vector'); 
                %Tri operation to ensure Hermitian symmetry of R_inv
                upper_tri_diag = triu(R_inv); %Obtain upper triangle
                upper_tri_pos_one = triu(R_inv,1); %Obtain upper triangle excluding diagonal elements 
                lower_tri_neg_one = upper_tri_pos_one'; %Reflect upper triangle to obtain lower triangle
                R_inv = upper_tri_diag + lower_tri_neg_one; %Add upper and lower triangle         

            end
            
            %Form reference phasor
            c = 1/(norm(rx_m_cap)^2); %Normalizing  factor
            V = c*(w_cur.')*input_cur; %Form V(k+1)
            w_prev1 = w_cur; %Form w(k-1)
            
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
%Compute the ensemble average of BER vs power per symbol
ave_BER = mean(BER,2); %Averaging along each row  

%Compute power value at specified BER value
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
save_loc = [pathname,filename,'_',num2str(M,'%u'),rp.format,'_BR_',num2str(R,'%10.3G'),'_FL_',num2str(N,'%u'),'_LLW_',num2str(laser_linewidth,'%10.3G'),'_FO_',num2str(frequency_offset,'%10.3G'),'_NA_',num2str(rp.NA,'%u'),'.mat'];
save(save_loc);

%%
%Get the screensize to specify figure size and location
scrsz = get(0,'ScreenSize'); 

% Specify position of figure on screen. rect = [left, bottom, width, height]
figure('OuterPosition',[1 scrsz(4)/2 scrsz(3)/4 scrsz(4)/2],'Name','BER vs Launch power') 
%Plot the BER vs power graph
semilogy(rp.power,ave_BER,'-b.','linewidth',1,'markerfacecolor','r')
title('BER vs Launch power'),
xlabel('Launch power (dBm)'),ylabel('BER'),grid;    

%%
%Find the elapsed time for the code to run
time2 = clock;
elapsed_minutes = etime(time2, time1)/60;
fprintf('Elapsed Time = %6.2f minutes\n', elapsed_minutes);
disp('End of program run');
disp('------------------');

end