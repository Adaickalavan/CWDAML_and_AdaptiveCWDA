function [rp] = NL_par(rp)

% Compute the nonlinear parameter list
h = 6.62606957e-34; %Planck constant, Units: J.s
C = 2.99792458e8; %Speed of light, Units: m/s

alpha = rp.alpha/4.343; %alpha = rp.alpha(dB/km)/4.343, Units: 1/km
% alpha = 10^(rp.alpha/10); (WRONG) %Convert units of alpha, Units: 1/km

L_eff = (1 - exp(-alpha*rp.L))/alpha; %Units: km
% L_eff = 1/alpha; 

rp.coeff = rp.gamma*L_eff; %Nonlinear coefficient, Units: rad/W 

nsp = 0.5*10^(rp.NF/10.0); %Spontaneous emission factor of optical amplifier, Units: arbitrary  
S_ASE = nsp*(10^(rp.G/10.0)-1)*h*C/(1e-9*rp.lambda0); %Spectral density of amplifier sponataneous emission noise, Units: J
rp.P_ASE = S_ASE*rp.Bo; %Noise power in single polarization within optical amplifier bandwidth/signal bandwidth, Units: W

end
