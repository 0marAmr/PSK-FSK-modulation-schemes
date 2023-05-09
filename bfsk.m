%{
*  FILE DESCRIPTION
*  -------------------------------------------------------------------------------------------------------------------
*  File:         Project3.m
*
*  Description:   Signal space representation and BER calculations for BFSK modulation schemes
*
*  -------------------------------------------------------------------------------------------------------------------
*  Author:       Omar Amr, Reem Mohamed
*  Date:         8 may 2023
%}
clc;
clear all;
close all;
%% Declarations 
no_of_bits = 1e5 + 1;
M = 10;              % number of realizations in the ensamble
N = 11;               % No of SNRs (channels) 

%%% BER declrations %%%
data_bits = randi([0 1],M,no_of_bits);  % data used to create the ensamble
data_sent = data_bits(1,:);             % bitstream to calculate BER
BFSK_demapped = zeros(N,no_of_bits);    % de-mapper output
Ch_out = zeros(N,no_of_bits);

%%% ensamble declrations %%%
no_of_samples = 10;
Tb = 1e-3;
Eb = 1; 

mapped_ensamble = zeros(M, no_of_samples * no_of_bits);
mapped_ensamble_delayed = [];

% Base-band signals
t = Tb/no_of_samples:(1/no_of_samples)*Tb:Tb*(1);
s1_bb = repelem(sqrt((2*Eb)/Tb),no_of_samples);                 % mapping zero
s2_bb = sqrt((2*Eb)/Tb) * exp(2*1i*pi*t/Tb);                    % mapping one


%% Mapper 
%%% input bitstream mapping %%%
FSK_mapped_data = data_sent;
for bit = 1 : no_of_bits
    if (FSK_mapped_data(bit) == 0)
        FSK_mapped_data(bit) = 1i;
    end
end

%%% Ensamble Mapping %%%
for realization = 1 : M
    for bit_count = 1 : no_of_bits
        if(data_bits(realization,bit_count) == 1)
            mapped_ensamble(realization,(bit_count-1)*no_of_samples + 1 : no_of_samples *bit_count) = s2_bb;
        else
             mapped_ensamble(realization,(bit_count-1)*no_of_samples + 1 : no_of_samples *bit_count) = s1_bb;
        end
    end
    
    %%% adding random delay to the ensamble %%%
    delay = randi([1,no_of_samples]) + 1; % delay = 1 means there is no delay, hence we add a +1 offset
    mapped_ensamble_delayed = vertcat(mapped_ensamble_delayed, mapped_ensamble(realization, delay : end - no_of_samples + (delay -1) ));
end

%% Channel
SNR_dB = 0:1:N;                             % different values for SNR in dB
SNR= 10.^(SNR_dB./10);                      % linear scale SNR = EB/No =1/No 
Eb = 1; 
No = Eb./SNR;

for CH_i = 1 : N                                
    I_noise =  randn(1,no_of_bits);
    Q_noise =  randn(1,no_of_bits);
    Ch_out(CH_i,:) =  (real(FSK_mapped_data) + sqrt(No(CH_i)/2) * I_noise)  + ...
         1i * (imag(FSK_mapped_data) + sqrt(No(CH_i)/2) * Q_noise);
end

%% demapper
for CH_i = 1: N
    for i = 1 : (no_of_bits)
        signal_arg = angle(Ch_out(CH_i,i));
        if(signal_arg<0) 
            signal_arg = signal_arg + 2*pi;
        end
        if(signal_arg <= pi/4 || signal_arg > 5*pi/4)
            BFSK_demapped(CH_i,i) = 1;
        else
            BFSK_demapped(CH_i,i) = 0;
        end
    end
end

%% BER Calculator
BFSK_pract_BER = zeros(1,N);
BFSK_Theo_BER = zeros(1,N);
for CH_i = 1: N
    % BFSK BER calculations
    [~,BFSK_pract_BER(CH_i)] = symerr(BFSK_demapped(CH_i,:),data_bits(1,:)); 
     BFSK_Theo_BER(CH_i) = 1/2.*erfc(sqrt( 1/(2*No(CH_i) ))); 
end

%% PSD calculation 
for column = 1: size(mapped_ensamble_delayed,2)
    R_tau(column) =(1/M) * sum(mapped_ensamble_delayed(:, 1) .* (mapped_ensamble_delayed(:, column)));
end
% R_tau = [fliplr(stat_autocorrelation(2:end)) stat_autocorrelation]; %Flipping R_tau to the -ve quad

%% Opututs

%%% BFSK constellation plot %%%
plot(real(Ch_out(N,:)),imag(Ch_out(N,:)),'gx', real(1), imag(0),'ko', real(0), imag(1i),'ko')
title('Signal-Space');
hold on;
x = -1:0.1:2;  % Define the range of x values from 0 to 2 with a step size of 0.01
y = x;  
plot(x,y,'--');

%%% BER Plot %%%
figure; 
EbNo_dB = 1:1:N;            % x-axis
EbNo = 10.^(EbNo_dB/10);    % y-axis
semilogy(EbNo_dB - 1,BFSK_Theo_BER,'--');   % BFSK BER theoretical Plot
grid on;
ylabel('BER');
xlabel('E_b/N_0');
title('BFSK modulation scheme BER plot')
hold on;
semilogy(EbNo_dB - 1,BFSK_pract_BER,'*-');  % BFSK BER experimental Plot
legend({'BFSK Theo BER','BFSK Exp BER'},'Location','southwest');

%%%% plottion PSD %%%

L = length(R_tau);    % Length of signal
n = 2^nextpow2(L);      
T = 0.1*Tb;           % Sampling period
Fs = 1/T;             % Sampling frequency                      
t = (0:L-1)*T;        % Time vector

%Compute the Fourier transform of the signal.
PSD = fft(R_tau,n);
PSD = fftshift(PSD);
f = Fs * (-n/2:n/2-1)/n;    % frequency (x-axis)

%Compute the two-sided spectrum P2
P2 = abs(PSD/n).^2;

Smoothed_PSD = smooth(P2); %Smooth the curve from the noise, to be plotted
figure;
subplot(1,1,1,'LineWidth',5)
plot(f,Smoothed_PSD) 
title("PSD")
xlabel("f(Hz)")
ylabel("S(f)")

