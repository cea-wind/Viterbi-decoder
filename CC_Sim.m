%% Convolution codes encode and decode simulation
%
%    Autor: Cao Wenhui
%    Last Modify:2015-07-23
%    Runtime:MATLAB(R) 2014a
%

%% 
% clear everything 
clc;clear all;close all;

%%
% Simulation parameter setting
EbN0_dB = 2:0.5:4.5;                   
FRAMES_NUM = 5000;                         
trel = poly2trellis(7, [171 133]);
%% BER Monte Carlo Simulation
bitError = zeros(size(EbN0_dB));
bitError2 = zeros(size(EbN0_dB));
RATE=1/2;
%%
% BER in different EbN0
for nEbN0 = 1:length(EbN0_dB)
%%
% Monte Carlo Simulation
    for nF=1:FRAMES_NUM
        %%
        % encode
        message = randi([0 1],1,100);
        encodeData = convenc(message,trel);
        %%
        % modulate
        transmitSignal = 1-2*encodeData ;   % 0-1;1--1

        %%
        % AWGN Channel, the relationship between SNR and EbN0
        
        SNR_dB = EbN0_dB((nEbN0)) + 10*log10(2)+10*log10(RATE);
        SNR = 10^(SNR_dB/10);
        noise = randn(1,length(transmitSignal));
        noise = noise/sqrt(SNR);     
        receiveSignal = transmitSignal + noise;

        %%
        % decode
        tblen = 30;
        decoded1 = vitbiDecoder(receiveSignal,trel,tblen);
        %% Comparing with vitdec,
        % The decode setting of vitdec didn't point out clearly, results of
        % two function may different.
        decoded2 = vitdec(receiveSignal,trel,tblen,'cont','unquant');
        if(sum(abs(decoded1-decoded2))~=0)
             sum(abs(message(1:end-tblen)-decoded1(tblen+1:end)))
             sum(abs(message(1:end-tblen)-decoded2(tblen+1:end)))
            disp('alert');
        end
        
        %%
        % BER,FER and iterations
        bitError(nEbN0) = bitError(nEbN0) + ...
            sum(abs(message(1:end-tblen)-decoded1(tblen+1:end)));
        bitError2(nEbN0) = bitError2(nEbN0) + ...
            sum(abs(message(1:end-tblen)-decoded2(tblen+1:end)));

    end
end
BER = bitError/5000/70;
BER2 = bitError2/5000/70;