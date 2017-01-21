# CWDAML and Adaptive CWDA

I had proposed and implemented 3 novel adaptive learning algorithms for time series data prediction and classification utilising modified linear and logistic regression techniques for optical communication receiver designs in my Ph.D. endeavour at National University of Singapore, Singapore, in the period Aug 2010 - Jul 2014. 

This repository contains the Matlab code to simulate [1] Complex-Weighted Decision-Aided Maximum-Likelihood (CWDAML) and [2] Adaptive Complex-weighted Decision-Aided (Adaptive CWDA), algorithms for M-ary Phase-Shift Keying (MPSK) and M-ary Quadrature Amplitude Modulation (MQAM) modulation formats in the presence of laser phase noise and frequency offset impairments at the receiver. The simulations consists of transmitter, laser phase noise, frequency offset between transmitter and receiver oscillator lasers, IQ receiver, baseband digital signal processing algorithm comprising CWDAML or adaptive CWDA algorithm.

For explanation of algorithms:
- Please refer to my LinkedIn profile for links to my published journal/conference papers where the algorithms have been explained in detail.
- LinkedIn: https://www.linkedin.com/in/adaickalavan

For CWDAML:
- The codes will draw BER vs SNR graph.
- In the Matlab folder, you only have to modify the parameters in “Batch_file_16QAM.m” or "Batch_file_4PSK.m" which is the main file. 
- The batch file will call all other files which are function files.
- For a start with my Matlab codes, you can simply run the batch file as-is without any modifications. It will produce a BER vs SNR graph corresponding to that of CWDAML algorithm. 
- The simulation will take several minutes to complete.
- You will need to create a "Simulation Results" folder in the same directory as the rest of the codes to store the results at the end of the simulation.
- The default settings in "Batch_file_4PSK.m" are: 
  - symbol rate = 28e9 symbols/s
  - Filter length = 15
  - Laser linewidth = 11.2e6 Hz
  - Frequency offset = 2.8e9 Hz
- The default settings in "Batch_file_16QAM.m" are: 
  - symbol rate = 28e9 symbols/s
  - Filter length = 12
  - Laser linewidth = 1.12e6 Hz
  - Frequency offset = 2.8e9 Hz  

For Adaptive CWDA:
- Codes will be added soon.
