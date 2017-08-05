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
rp.No_sent_symbols = 2e6; %Total no of symbols to send before computing BER
rp.total_error_bit = 500; %Total no bit errors to count before computing BER

%Fuction parameters
rp.idl_dcs_fdb = 0; %Ideal decision feedback 

%Call the functions     
rp.M = 4;
rp.format = 'PSK'; 
rp.bit_rate = log2(rp.M)*28e9; %Bit rate, Units: bits/s

%Non-linear parameter list
C = 2.99792458e8; %Speed of light, Units: m/s
rp.lambda0 = 1550; %Wavelength, Units: nm
% Bo_wavelength = 0.5; %optical filter bandwidth in wavelength, Units: nm
% rp.Bo = Bo_wavelength*C*1e9/(rp.lambda0^2); %Optical amplifier bandwidth, Units: Hz
rp.Bo = rp.bit_rate/log2(rp.M); %Optical amplifier bandwidth, Units: Hz
rp.gamma = 1.2; %Fiber nonlinearity coefficient, Units: rad/(W.km)
rp.alpha = 0.2; %Fiber loss, Units: dB/km
rp.L = 100; %Fiber span length, Units: km
rp.G = rp.alpha*rp.L; %Amplifier gain, Units: dB
rp.NF = 4.5; %Amplifier noise figure such that n_sp = 1.41, Units: dB

%Other parameters
LSR = 200e3/28e9; %200kHz laser linewidth per symbol;
FSR = 0; %Frequency/Symbol_rate

LBR = (1/log2(rp.M))*LSR; %Linewidth/Bit_rate
LLW = rp.bit_rate*LBR;
rp.laser_linewidth = LLW;
FBR = (1/log2(rp.M))*FSR; %Frequency/Bit_rate
FO = rp.bit_rate*FBR; %Frequency offset
rp.frequency_offset = FO; %Frequency offset, Units: Hz

LP = (-10:2:10).';
rp.read_BEP = 2.5E-2;
rp.NA = 41;
rp.power = LP; %signal power 

% rp.training_length = 300; %Training_length. Not applicable to Block_Mth_power scheme.  
% rp.freq_est_length = 4e3; %No_sent_symbols to be used for frequency estimation   
% value = NL_Adaptive_CWDA_MPSK(rp);
 
rp.filter_length = 21;
rp.training_length = 100; %Training_length. Not applicable to Block_Mth_power scheme.  
rp.freq_est_length = 4e3; %No_sent_symbols to be used for frequency estimation    
value = NL_CW_recursiveW_DAML_MPSK(rp);

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