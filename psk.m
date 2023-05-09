%{
*  FILE DESCRIPTION
*  -------------------------------------------------------------------------------------------------------------------
*  File:         Project3.m
*
*  Description:   Signal space representation and BER calculations for PSK modulation schemes
*
*  -------------------------------------------------------------------------------------------------------------------
*  Author:       Omar Amr, Reem Mohamed
*  Date:         8 may 2023
%}
clc;
clear all;
%% Declarations 
no_of_bits = 1.8e5;
data_bits = randi([0 1],1,no_of_bits);  % i/p stream of bits (row vector)
N = 11;      % No of SNRs (channels) 

%%% Channel output declrations %%%
BPSK_Ch_out = zeros(N,no_of_bits); %group of one bit-->BPSK
QPSK_Ch_out = zeros(N,no_of_bits/2); %group of 2 bits-->QPSK
QPSK_bin_Ch_out = zeros(N,no_of_bits/2); %group of 2 bits-->QPSK (binary encoded)
M8PSK_Ch_out = zeros(N,no_of_bits/3); %group of 3 bits-->M_8PSK
QAM16_Ch_out = zeros(N,no_of_bits/4); %group of 4 bits-->QAM16

%%% Channel demapper declrations %%%
BPSK_demapped = zeros(N,no_of_bits);    %demapping to stream of single bit to compare with original data bits and calculate BER practically
QPSK_demapped = zeros(N,no_of_bits); %demapping to stream of single bit tocompare with original data bits and calculate BER practically
QPSK_bin_demapped = zeros(N,no_of_bits); %demapping to stream of single bit tocompare with original data bits and calculate BER practically
M8PSK_demapped = zeros(N,no_of_bits); %demapping to stream of single bit to compare with original data bits andcalculate BER practically
QAM16_demapped = zeros(N,no_of_bits); %demapping to stream of single bit to compare with original data bits and calculate BER practically

%% Mapper
% BPSK mapping
BPSK_mapped_data = 2 *data_bits - 1;

% QPSK mapping
% A mapper that converts groups of two bits into grey code
QPSK_mapper = [-1-1i, -1+1i, 1-1i, 1+1i];          
QPSK_binary_reshaped= reshape(data_bits,2,[])';    % data reshaped into two rows (QPSK Grouping)
QPSK_mapped_data = QPSK_mapper(bi2de(QPSK_binary_reshaped,'left-msb')+1);

% QPSK_bin mapping
% A mapper that converts groups of two bits not into grey code
QPSK_mapper_bin = [-1-1i, -1+1i, 1+1i, 1-1i];          
QPSK_binary_reshaped_bin= reshape(data_bits,2,[])';    % data reshaped into two rows (QPSK Grouping)
QPSK_bin_mapped_data = QPSK_mapper_bin(bi2de(QPSK_binary_reshaped_bin,'left-msb')+1);

% 8PSK mapping
% A mapper that converts groups of two bits into grey code
M8PSK_mapper = [1, exp(1i*(pi/4)), exp(1i*(3*pi/4)), 1i, exp(1i*(7*pi/4)), -1i, -1, exp(1i*(5*pi/4))];      
M8PSK_binary_reshaped= reshape(data_bits,3,[])';    % data reshaped into three rows (8PSK Grouping)
M8PSK_mapped_data = M8PSK_mapper(bi2de(M8PSK_binary_reshaped,'left-msb')+1);

% 8PSK mapping
QAM16_mapper = [-3-3i, -3-1i, -3+3i, -3+1i, -1-3i, -1-1i, -1+3i, -1+1i,... 
              3-3i,  3-1i,  3+3i,  3+1i, 1-3i,  1-1i,  1+3i,  1+1i];
QAM16_binary_reshaped= reshape(data_bits,4,[])';    % data reshaped into four rows (QAM16 Grouping)
QAM16_mapped_data = QAM16_mapper(bi2de(QAM16_binary_reshaped,'left-msb')+1);

%% Channel
% average energy per bit
Eb_BPSK = 1; %costellation point is root(EB) =1 --> EB=1                             
Eb_QPSK = ((4*2)/4)/2;  %or costellation point is root(E/2) =1 --> E=2 -->E=2EB --> EB=1 -->E energy of signal
Eb_QPSK_bin = ((4*2)/4)/2;  %or costellation point is root(E/2) =1 --> E=2 -->E=2EB --> EB=1 -->E energy of signal
Eb_M8PSK = ((8*1)/8)/3; %or costellation point on Q1 is root(E) =1 --> E=1 -->E=3EB --> EB=1/3>E energy of signal
Eb_QAM16 =((1/16*4*(2+10+10+18))/4); %avg energy per symbol and symbol is 4 bits -->E=4EB

SNR_dB = 0:1:N;                             % different values for SNR in dB
SNR= 10.^(SNR_dB./10);                      % linear scale SNR = EB/No =1/No 

% Noise coff. calculations
No_BPSK = Eb_BPSK./SNR;  
No_QPSK = Eb_QPSK./SNR;
No_QPSK_bin = Eb_QPSK_bin./SNR;
No_M8PSK = Eb_M8PSK./SNR; 
No_QAM16 = Eb_QAM16./SNR; 

for i = 1 : N                            
    % channel effect on BPSK scheme
    BPSK_Ch_out(i,:) = BPSK_mapped_data + sqrt(No_BPSK(i)/2) * randn(1,no_of_bits); 
    
    % channel effect on QPSK scheme
    I_noise_QPSK =  randn(1,no_of_bits/2);
    Q_noise_QPSK =  randn(1,no_of_bits/2);
    QPSK_Ch_out(i,:) =  (real(QPSK_mapped_data) + sqrt(No_QPSK(i)/2) * I_noise_QPSK)  + ...
         1i * (imag(QPSK_mapped_data) + sqrt(No_QPSK(i)/2) * Q_noise_QPSK);

    % channel effect on QPSK_NG scheme
    I_noise_QPSK_bin =  randn(1,no_of_bits/2);
    Q_noise_QPSK_bin =  randn(1,no_of_bits/2);
    QPSK_bin_Ch_out(i,:) =  (real(QPSK_bin_mapped_data) + sqrt(No_QPSK_bin(i)/2) * I_noise_QPSK_bin)  + ...
         1i * (imag(QPSK_bin_mapped_data) + sqrt(No_QPSK_bin(i)/2) * Q_noise_QPSK_bin);
     
    % channel effect on 8PSK scheme
    I_noise_M8PSK =  randn(1,no_of_bits/3);
    Q_noise_M8PSK =  randn(1,no_of_bits/3);
    M8PSK_Ch_out(i,:) =  (real(M8PSK_mapped_data) + sqrt(No_M8PSK(i)/2) * I_noise_M8PSK)  + ...
         1i * (imag(M8PSK_mapped_data) + sqrt(No_M8PSK(i)/2) * Q_noise_M8PSK);
     
    % channel effect on QAM16 scheme
    I_noise_QAM16 =  randn(1,no_of_bits/4);
    Q_noise_QAM16 =  randn(1,no_of_bits/4);
    QAM16_Ch_out(i,:) =  (real(QAM16_mapped_data) + sqrt(No_QAM16(i)/2) * I_noise_QAM16)  + ...
         1i * (imag(QAM16_mapped_data) + sqrt(No_QAM16(i)/2) * Q_noise_QAM16);
end

%% Demapper and BER Calculator
QPSK_demapped_symbols = zeros(no_of_bits/2,2);
QPSK_bin_demapped_symbols= zeros(no_of_bits/2,2);
M8PSK_demapped_symbols = zeros(no_of_bits/3,3);
QAM16_demapped_symbols = zeros(no_of_bits/4,4);

%%% declrations for 8PSK %%%
M8PSK_decision_angles = pi/8:pi/4:((15*pi)/8) ;   % decision angles midway constellation points
% generate 3-bit grey code to map the received messages to in case of8PSK
decimal_grey_code = [0 1 3 2 6 7 5 4];
grey_code_3b = de2bi([0 1 3 2 6 7 5 4],'left-msb');

for CH_i = 1: N
    % BPSK de-mapper
    for i = 1 : no_of_bits
        if (BPSK_Ch_out(CH_i,i) > 0)
            BPSK_demapped(CH_i,i) = 1;
        else
            BPSK_demapped(CH_i,i) = 0;
        end
    end

    % QPSK de-mapper
    for i = 1 : (no_of_bits/2)
        if (real(QPSK_Ch_out(CH_i,i)) > 0 && imag(QPSK_Ch_out(CH_i,i)) > 0 )
            QPSK_demapped_symbols(i,:) = [1 1];
        elseif (real(QPSK_Ch_out(CH_i,i)) > 0 && imag(QPSK_Ch_out(CH_i,i)) < 0 )
            QPSK_demapped_symbols(i,:) = [1 0];
        elseif (real(QPSK_Ch_out(CH_i,i)) < 0 && imag(QPSK_Ch_out(CH_i,i)) < 0 )
            QPSK_demapped_symbols(i,:) = [0 0];
        else
            QPSK_demapped_symbols(i,:) = [0 1];
        end
    end
    % reshaping the demapped data into output bit stream (a row vector)
    QPSK_demapped(CH_i,:) = reshape(QPSK_demapped_symbols',1,[]);
    
    % QPSK_bin de-mapper
    for i = 1 : (no_of_bits/2)
        if (real(QPSK_bin_Ch_out(CH_i,i)) > 0 && imag(QPSK_bin_Ch_out(CH_i,i)) > 0 )
            QPSK_bin_demapped_symbols(i,:) = [1 0];
        elseif (real(QPSK_bin_Ch_out(CH_i,i)) > 0 && imag(QPSK_bin_Ch_out(CH_i,i)) < 0 )
            QPSK_bin_demapped_symbols(i,:) = [1 1];
        elseif (real(QPSK_bin_Ch_out(CH_i,i)) < 0 && imag(QPSK_bin_Ch_out(CH_i,i)) < 0 )
            QPSK_bin_demapped_symbols(i,:) = [0 0];
        else
            QPSK_bin_demapped_symbols(i,:) = [0 1];
        end
    end
    % reshaping the demapped data into output bit stream (a row vector)
    QPSK_bin_demapped(CH_i,:) = reshape(QPSK_bin_demapped_symbols',1,[]);
    
    % 8PSK de-mapper
    % Algoritm idea by: Eng. Reem M. Ashour =) <3
    for i = 1 : (no_of_bits/3)
        % by default, we assume that the msg belongs to constellation 000
        M8PSK_demapped_symbols(i,:) = grey_code_3b(1,:);
        signal_arg = angle(M8PSK_Ch_out(CH_i,i));
        if(signal_arg<0) 
            signal_arg = signal_arg + 2*pi;
        end
        for arg = 2 : length(M8PSK_decision_angles) 
            % check for the argument of the complex signal sent
            if (signal_arg >= M8PSK_decision_angles(arg - 1) && signal_arg <= M8PSK_decision_angles(arg))
                M8PSK_demapped_symbols(i,:) = grey_code_3b(arg,:);
                break;
            end
        end
    end
    M8PSK_demapped(CH_i,:) = reshape(M8PSK_demapped_symbols',1,[]);

    % QAM16 de-mapper
    for i = 1 : (no_of_bits/4)
        if (real(QAM16_Ch_out(CH_i,i)) > 0 && real(QAM16_Ch_out(CH_i,i)) < 2 && imag(QAM16_Ch_out(CH_i,i)) > 0 && imag(QAM16_Ch_out(CH_i,i)) <2 )
            QAM16_demapped_symbols(i,:) = [1 1 1 1];
        elseif (real(QAM16_Ch_out(CH_i,i)) > 0 && real(QAM16_Ch_out(CH_i,i)) < 2 && imag(QAM16_Ch_out(CH_i,i)) <0 && imag(QAM16_Ch_out(CH_i,i)) >-2 )
            QAM16_demapped_symbols(i,:) = [1 1 0 1];
        elseif (real(QAM16_Ch_out(CH_i,i)) > 0 && real(QAM16_Ch_out(CH_i,i)) < 2 && imag(QAM16_Ch_out(CH_i,i)) >2 )
            QAM16_demapped_symbols(i,:) = [1 1 1 0];
        elseif (real(QAM16_Ch_out(CH_i,i)) > 0 && real(QAM16_Ch_out(CH_i,i)) < 2 && imag(QAM16_Ch_out(CH_i,i)) <-2 )
            QAM16_demapped_symbols(i,:) = [1 1 0 0];
        elseif (real(QAM16_Ch_out(CH_i,i)) > 2 && imag(QAM16_Ch_out(CH_i,i)) > 0 && imag(QAM16_Ch_out(CH_i,i)) <2 )
            QAM16_demapped_symbols(i,:) = [1 0 1 1];
        elseif (real(QAM16_Ch_out(CH_i,i)) > 2 && imag(QAM16_Ch_out(CH_i,i)) < 0 && imag(QAM16_Ch_out(CH_i,i)) >-2 )
            QAM16_demapped_symbols(i,:) = [1 0 0 1];
        elseif (real(QAM16_Ch_out(CH_i,i)) > 2 && imag(QAM16_Ch_out(CH_i,i)) > 2 )
            QAM16_demapped_symbols(i,:) = [1 0 1 0];
        elseif (real(QAM16_Ch_out(CH_i,i)) > 2 && imag(QAM16_Ch_out(CH_i,i)) <- 2 )
            QAM16_demapped_symbols(i,:) = [1 0 0 0];    
        elseif (real(QAM16_Ch_out(CH_i,i)) < 0 && real(QAM16_Ch_out(CH_i,i)) >-2 && imag(QAM16_Ch_out(CH_i,i)) > 0 && imag(QAM16_Ch_out(CH_i,i)) <2 )
            QAM16_demapped_symbols(i,:) = [0 1 1 1];
        elseif (real(QAM16_Ch_out(CH_i,i)) < 0 && real(QAM16_Ch_out(CH_i,i)) >-2 && imag(QAM16_Ch_out(CH_i,i)) < 0 && imag(QAM16_Ch_out(CH_i,i)) >-2 )
            QAM16_demapped_symbols(i,:) = [0 1 0 1];
        elseif (real(QAM16_Ch_out(CH_i,i)) < 0 && real(QAM16_Ch_out(CH_i,i)) >-2 && imag(QAM16_Ch_out(CH_i,i)) >2 )
            QAM16_demapped_symbols(i,:) = [0 1 1 0];
        elseif (real(QAM16_Ch_out(CH_i,i)) < 0 && real(QAM16_Ch_out(CH_i,i)) >-2 && imag(QAM16_Ch_out(CH_i,i)) <-2 )
            QAM16_demapped_symbols(i,:) = [0 1 0 0];
        elseif (real(QAM16_Ch_out(CH_i,i)) < -2 && imag(QAM16_Ch_out(CH_i,i)) > 0 && imag(QAM16_Ch_out(CH_i,i)) <2 )
            QAM16_demapped_symbols(i,:) = [0 0 1 1];
        elseif (real(QAM16_Ch_out(CH_i,i)) < -2 && imag(QAM16_Ch_out(CH_i,i)) < 0 && imag(QAM16_Ch_out(CH_i,i)) >-2 )
            QAM16_demapped_symbols(i,:) = [0 0 0 1];
        elseif (real(QAM16_Ch_out(CH_i,i)) < -2 && imag(QAM16_Ch_out(CH_i,i)) >2 )
            QAM16_demapped_symbols(i,:) = [0 0 1 0];
        else
            QAM16_demapped_symbols(i,:) = [0 0 0 0];
        end
    end
     QAM16_demapped(CH_i,:) = reshape(QAM16_demapped_symbols',1,[]);%to get stream of only single bit to compare with original bits and get BER practically
    
    
end

%% BER Calculator
BPSK_pract_BER = zeros(1,N);
QPSK_pract_BER = zeros(1,N);
QPSK_bin_pract_BER = zeros(1,N);
M8PSK_pract_BER = zeros(1,N);
QAM16_pract_BER = zeros(1,N);

BPSK_Theo_BER = zeros(1,N);
QPSK_Theo_BER = zeros(1,N);
M8PSK_Theo_BER = zeros(1,N);
QAM16_Theo_BER = zeros(1,N);

for CH_i = 1: N
    % BPSK calculations
    [~,BPSK_pract_BER(CH_i)] = symerr(BPSK_demapped(CH_i,:),data_bits); %compare between mapped and actual bits and calculate BER   
    BPSK_Theo_BER(CH_i) = 1/2.*erfc(sqrt(1/No_BPSK(CH_i)));%EB=1
    
    % QPSK calculations
    [~,QPSK_pract_BER(CH_i)] = symerr(QPSK_demapped(CH_i,:),data_bits); %compare between mapped and actual bits and calculate BER    
    QPSK_Theo_BER(CH_i) = 1/2.*erfc(sqrt(1/No_QPSK(CH_i)));%SER =erf(EB/No)-->BER=SER/2-->EB=1
    
    % QPSK_bin calculations
    [~,QPSK_bin_pract_BER(CH_i)] = symerr(QPSK_bin_demapped(CH_i,:),data_bits); %compare between mapped and actual bits and calculate BER  
    
    % 8PSK calculations
    [~,M8PSK_pract_BER(CH_i)] = symerr(M8PSK_demapped(CH_i,:),data_bits); %compare between mapped and actual bits and calculate BER   
    M8PSK_Theo_BER(CH_i) = erfc(sqrt(1/No_M8PSK(CH_i))*sin(pi/8))/3;%SER =erf(root(E/No)sin(pi/8))-->E=1-->BER=SER/3-->grey coding
     
    % QAM16 calculations
    [~,QAM16_pract_BER(CH_i)] = symerr(QAM16_demapped(CH_i,:),data_bits); %compare between mapped and actual bits and calculate BER  
     QAM16_Theo_BER(CH_i) = 3/2.*erfc(sqrt(1/No_QAM16(CH_i)))/4; %SER=3/2*erf(root(EO/NO))-->EO=1-->BER=SER/4-->grey coding
end


%% outputs
close all;
%%%% Signal-space plot %%%
% BPSK constellation 
subplot(2,2,1,'LineWidth',2)
plot(real(BPSK_Ch_out(N,:)),imag(BPSK_Ch_out),'rx', real(1), imag(0),'ko', real(-1), imag(0),'ko')
title('BPSK Signal-Space');
% QPSK constellation 
subplot(2,2,2,'LineWidth',2)
plot(real(QPSK_Ch_out(N,:)),imag(QPSK_Ch_out(N,:)),'gx', real(QPSK_mapper), imag(QPSK_mapper),'ko')
title('QPSK Signal-Space');
% 8PSK constellation 
subplot(2,2,3,'LineWidth',2)
plot(real(M8PSK_Ch_out(N,:)),imag(M8PSK_Ch_out(N,:)),'yx', real(M8PSK_mapper), imag(M8PSK_mapper),'ko')
title('8PSK Signal-Space');
% QAM16 constellation 
subplot(2,2,4,'LineWidth',2)
plot(real(QAM16_Ch_out(N,:)),imag(QAM16_Ch_out(N,:)),'bx', real(QAM16_mapper), imag(QAM16_mapper),'ko')
title('QAM16 Signal-Space');


%%% BER plot for all modulation schemes %%%
figure;
EbNo_dB = 1:1:N;            % x-axis
EbNo = 10.^(EbNo_dB/10);    % y-axis

% BPSK Plot
semilogy(EbNo_dB - 1,BPSK_Theo_BER,'--','Color',[1 0 0]);
grid on;
ylabel('BER');
xlabel('E_b/N_0');
title('4 Modulation Schemes BER')
hold on;
semilogy(EbNo_dB - 1,BPSK_pract_BER,'*-','Color',[1 0 0]);
hold on; 
% QPSK Plot
semilogy(EbNo_dB - 1,QPSK_Theo_BER,'--','Color',[0 1 0]);
hold on; 
semilogy(EbNo_dB - 1,QPSK_pract_BER,'x-','Color',[0 1 0]);
hold on; 
% 8PSK Plot
semilogy(EbNo_dB - 1,M8PSK_Theo_BER,'--','Color',[0 0 1]);
hold on; 
semilogy(EbNo_dB - 1,M8PSK_pract_BER,'.-','Color',[0 0 1]);
hold on; 
% QAM16 Plot
semilogy(EbNo_dB - 1,QAM16_Theo_BER,'--','Color',[1 0 1]);
hold on; 
semilogy(EbNo_dB - 1,QAM16_pract_BER,'o-','Color',[1 0 1]);
legend({'BPSK Theo BER','BPSK Exp BER', 'QPSK Theo BER', 'QPSK Exp BER', '8PSK Theo BER',...
    '8PSK Exp BER','QAM16 Theo BER', 'QAM16 Exp BER'},'Location','southwest','NumColumns',4);

%%% comparing BER for grey and binary encoded QPSK %%%
figure; 
semilogy(EbNo_dB - 1,QPSK_Theo_BER,'--','Color',[0.5 0.5 0.5]); % QPSK theoretical 
title('grey vs binary encoded QPSK BER')
ylabel('BER');
xlabel('E_b/N_0');
grid on;
hold on; 
semilogy(EbNo_dB - 1,QPSK_pract_BER,'Color',[1 0 0]);           % QPSK grey encoded practical 
hold on; 
semilogy(EbNo_dB - 1,QPSK_bin_pract_BER,'x-','Color',[0 0 1]);      % QPSK binary encoded practical
legend({'QPSK Theo', 'QPSK grey', 'QPSK bin'},'Location','southwest','NumColumns',3);

% figure;
% subplot(4,1,1,'LineWidth',2)
% stairs(data_bits,'B','LineWidth',2);
% axis([1 100 -3 3])
% title('data bits');
% xlabel('time');
% ylabel('data');

% subplot(4,1,2,'LineWidth',2)
% stairs(BPSK_mapped_data,'B','LineWidth',2);
% axis([1 500 -3 3])
% title('mapped data');
% xlabel('time');
% ylabel('data');
% 
% subplot(4,1,3,'LineWidth',2)
% plot(BPSK_Ch_out(1),'B','LineWidth',2);
% axis([1 500 -3 3])
% title('channel out');
% xlabel('time');
% ylabel('data');
% 
% subplot(4,1,4,'LineWidth',2)
% stairs(BPSK_demapped,'B','LineWidth',2);
% axis([1 100 -3 3])
% title('demapped data bits');
% xlabel('time');
% ylabel('data');
% 


