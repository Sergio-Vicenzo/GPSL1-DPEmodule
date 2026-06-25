%----------------------------------------
%
% Acquisition plus fine acquisition
%
% HAYDEN 04.16.2025
%----------------------------------------

function acqResults = acquisition_hd (longsignal,settings,fineacqdata)

%% Initialization
% Carrier frequencies of detected signals
acqResults.carrFreq     = zeros(1, 32);
% C/A code phases of detected signals
acqResults.codePhase    = zeros(1, 32);
% Correlation peak ratios of the detected signals
acqResults.peakMetric   = zeros(1, 32); % here store snr
%acqResults.snr=zeros(1,32);
acqResults.doppler = zeros(1,32);

settings.Sample =  ceil(settings.samplingFreq*20);
sampleindex = 1: ceil(settings.samplingFreq*20);

%% Acquisition parameters
acq.prnList     = 1:32;	% PRN list
acq.freqStep    = 500;	% unit: Hz
acq.freqMin     = -10000;   % Minimum Doppler frequency
acq.freqNum     = 2*abs(acq.freqMin)/acq.freqStep+1;    % number of frequency bins
acq.datalen     = 20; % msec
acq.L           = 10; % number of ms to perform FFT

for freqband = 1 : acq.freqNum
    dopplershift = acq.freqMin + acq.freqStep*(freqband-1);
    carrier(freqband,:) = exp(1i*2*pi*(settings.IF+ dopplershift) * sampleindex ./ settings.samplingFreq);
end


fprintf('Acquiring GPS (');

for PRN = 1:32
    ocode = generateCAcode(PRN);%generate oringinal C/A code
    ocode = [ocode ocode];%in case oversize
    scode = ocode(ceil(sampleindex.*(settings.codeFreqBasis/settings.samplingFreq)));%sample C/A code
    correlation = zeros(acq.freqNum,settings.Sample);
    for idx = 1 : 20  % 20ms data
        for freqband = 1 : acq.freqNum %parallel code phase
            replica = scode;
            temp1 = longsignal(1+(idx-1)*settings.Sample:idx*settings.Sample).* carrier(freqband,:);
            temp2 = conj(fft(temp1));
            temp3 = fft(replica);
            correlation(freqband,:) = correlation(freqband,:) + abs(ifft(temp3.*temp2)).^2;
        end
    end
    [peak, fbin] = max(max(abs(correlation')));
    [peak, codePhase] = max(max(abs(correlation)));
    Doppler = acq.freqMin + acq.freqStep * (fbin-1);
    
    codechipshift = ceil(settings.samplingFreq/settings.codeFreqBasis);
    SNR = 10 * log10(peak^2/(sum(correlation(fbin,[1:codePhase-codechipshift codePhase+codechipshift:end]).^2)... % outside the main correlation peak
        /length(correlation(fbin,[1:codePhase-codechipshift codePhase+codechipshift:end]))));
    
    acqResults.peakMetric(PRN)=SNR; 
    if SNR >= 10  % acquisition thredhold
        fprintf('%02d ', PRN);
        acqResults.doppler(PRN)=Doppler; %raw doppler
        acqResults.codePhase(PRN)=codePhase-1;   
    else
         fprintf('. ');
    end
end
fprintf(')\n');



fprintf('Fine Acquiring GPS (');
%% Fine acquisition caculation, length = 10 ms
for PRN = 1 : 32
    if (acqResults.doppler(PRN))~=0
         fprintf('%02d ', PRN);
    else 
         fprintf('. ');
         continue;
    end
    caCode = generateCAcode(PRN);
    codeValueIndex = floor((1/settings.samplingFreq*(1:10*settings.Sample))/(1/settings.codeFreqBasis));
    longCaCode = caCode((rem(codeValueIndex, settings.codeLength) + 1));
    CarrSignal = fineacqdata((settings.Sample-acqResults.codePhase(PRN)):(settings.Sample-acqResults.codePhase(PRN))+10*settings.Sample - 1).* longCaCode;
    
    fftlength = length(CarrSignal) * 20;
    
%         fftSignal = abs(fft(CarrSignal, fftlength)); % For NSL STEREO L1 Band only
    fftSignal = abs(fftshift(fft(CarrSignal, fftlength)));  % For NSL STEREO LBand only

    halffftlength = ceil((fftlength)/2);
    [~, FreqPeakIndex] = max(fftSignal(1:halffftlength*2));
    fineDoppler = FreqPeakIndex * (settings.samplingFreq/fftlength);
    if settings.fileType == 2 % default to be 2
        %         if Acquired.Doppler(svindex)>= 0
        %             fineDoppler = -FreqPeakIndex * (signal.Fs/fftlength) + signal.Fs/2;
        %         else
        fineDoppler = -FreqPeakIndex * (settings.samplingFreq/fftlength) + settings.samplingFreq/2;
        %         end
        
        %         if Acquired.Doppler(svindex)>= 0
        %             fineDoppler = FreqPeakIndex * (signal.Fs/fftlength);
        %         else
        %             fineDoppler = -FreqPeakIndex * (signal.Fs/fftlength);
        %         end
    end
    %Acquired.fineFreq = [Acquired.fineFreq  fineDoppler];
    acqResults.carrFreq(PRN)=fineDoppler; %-settings.IF;
    % fprintf('SV[%2d] SNR = %2.2f, Code phase = %5d, Raw Doppler = %5d, Fine Doppler = %5f \n ', ...
    %     Acquired.sv(svindex),Acquired.SNR(svindex),Acquired.codedelay(svindex), ...
    %     Acquired.Doppler(svindex),Acquired.fineFreq(svindex)-signal.IF);
end
fprintf(')\n');
end