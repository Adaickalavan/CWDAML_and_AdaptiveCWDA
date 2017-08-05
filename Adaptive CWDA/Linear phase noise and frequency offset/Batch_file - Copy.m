%Karpaga Vinayagar - Pillaiyarpatti
%Adaickalavan Meiyappan - NUS

%%
%This is a standard batch file to execute functions
clear all
close all

disp('--------------------');
disp('Start of batch file');
time1 = clock; %Retrieve starting time to compute code run time later
format shortG %Number display format
save_loc = './Results/Batch Results.mat'; % Location to save batch file values

%%

%Initialize general run parameters (rp)   
stream0 = RandStream('mt19937ar','Seed',0);
rp.stream = stream0; %Random number stream
rp.total_run = 1; %Number of times to repeat simulation to obtain ensemble average
rp.BEP = 5e-4; %Bit error probability threshold until which to run the simulation
rp.read_BEP = 1e-3; %Read SNR at this bit error probability value. Note that: rp.read_BEP > rp.BEP
rp.SNR_begin = 6; %Signal to noise ratio per bit, Units: dB
rp.SNR_max = 20; %Maximum SNR limit before declaring a phase tracking failure, Units: dB
rp.step_size = 1; %Incremental value of SNR 0.5 dB or 1 dB, Units: dB
rp.No_sent_symbols = 5e6; %Total no of symbols to send before computing BER
rp.total_error_bit = 500; %Total no bit errors to count before computing BER
rp.idl_dcs_fdb = 0; %Ideal decision feedback
rp.quantize = 0; %Number of quantization bits, 0 = no_quantization

%Call the functions     
rp.M = 4;
rp.format = 'PSK';
symbolRate = 28e9;
LSR = 1.12e6/symbolRate; %Laser linewidth(Hz)/Symbol_Rate
FSR = 2.8e9/symbolRate; %Frequency offset(Hz)/Symbol_rate           

rp.bit_rate = log2(rp.M)*symbolRate; %Bit rate, Units: bits/s 
LBR = (1/log2(rp.M))*LSR; %Linewidth/Bit_rate
FBR = (1/log2(rp.M))*FSR; %Frequency/Bit_rate
FO = rp.bit_rate*FBR; %Frequency offset
LLW = rp.bit_rate*LBR; %Laser Linewidth
rp.laser_linewidth = LLW; %Laser linewidth, Units: Hz
rp.frequency_offset = FO; %Frequency offset, Units: Hz          

rp.training_length = 300; %Training_length. Not applicable to Block_Mth_power scheme.   
rp.freq_est_length = inf; %No_sent_symbols to be used for frequency estimation
value = Adaptive_CWDA_MPSK(rp);            

%%
%Save all variables from current workspace for later access
save(save_loc);

%%
%Find the elapsed time for the code to run
time2 = clock;
elapsed_minutes = etime(time2, time1)/60;
fprintf('Elapsed Time = %6.2f minutes\n', elapsed_minutes);
disp('End of batch file');
disp('------------------');
