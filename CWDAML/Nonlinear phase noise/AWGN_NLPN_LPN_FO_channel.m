function [rx, phase_noise, act_tpn] = AWGN_NLPN_LPN_FO_channel(m,phase_noise,rp)

%Received message
rx_NL = m;
block_length = length(rx_NL);

% %Nonlinear phase noise and AWGN impairment
% for i = 1:rp.NA
%     rx_NL = rx_NL + sqrt(rp.P_ASE/2)*(randn(block_length,1) + 1j*randn(block_length,1));
%     %Var['randn'] = 1 and Var[sqrt(rp.P_ASE/2)*(randn + 1j*randn)] = rp.P_ASE
%     rx_NL = rx_NL.*exp(1j*rp.coeff*abs(rx_NL).^2);
% end
% 
% %Linear phase noise and frequency offset impairment
% phase_noise = [sum(phase_noise) + rp.sigma*randn + rp.omega; rp.sigma*randn(block_length-1,1) + rp.omega];  %Current cumulative phase noise, Units: radians
% % symbol_no = (rp.noss + (0:block_length-1)).';
% % total_phase_noise = cumsum(phase_noise)+(rp.App/rp.freq_mod)*sin(rp.omega_FM*symbol_no);
% rx = rx_NL.*exp(1j*cumsum(phase_noise)); %Received signal with linear & nonlinear phase noise & AWGN                                            

%Actual accumulated nonlinear phase noise model given in paper
rx = zeros(block_length,1);
act_tpn = zeros(block_length,1);
for i = 1:block_length
    vector_1 = [rx_NL(i); sqrt(rp.P_ASE/2)*(randn(rp.NA,1) + 1j*randn(rp.NA,1))];
    vector_2 = abs(cumsum(vector_1)).^2;
    phi_NL = rp.coeff*sum(vector_2(2:end));
    phase_noise = phase_noise + rp.sigma*randn + rp.omega;
    E_NA = sum(vector_1);
    rx(i,1) = E_NA*exp(1j*phi_NL)*exp(1j*phase_noise);
    %Actual total phase noise
    act_tpn(i,1) = phi_NL + phase_noise;
end

end
