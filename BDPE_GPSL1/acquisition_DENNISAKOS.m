function acqResults = acquisition_DENNISAKOS(signal, settings)
%Function performs cold start acquisition on the collected "data". It
%searches for GPS signals of all satellites, which are listed in field
%"acqSatelliteList" in the settings structure. Function saves code phase
%and frequency of the detected signals in the "acqResults" structure.
%
%acqResults = acquisition(signal, settings)
%
%   Inputs:
%       signal        - raw signal from the front-end
%       settings      - Receiver settings. Provides information about
%                       sampling and intermediate frequencies and other
%                       parameters including the list of the satellites to
%                       be acquired.
%   Outputs:
%       acqResults    - Function saves code phases and frequencies of the
%                       detected signals in the "acqResults" structure. The
%                       field "carrFreq" is set to 0 if the signal is not
%                       detected for the given PRN number.

%--------------------------------------------------------------------------
%                           SoftGNSS v3.0
%
% Copyright (C) Dennis M. Akos
% Written by Sirish Jetti
% Updated by Xiaofan Li
% Should be noted that there is a potential code phase samples deviation in
% the results of acquisition, the largest deviation can be up to 8 samples
% or 0.5 chip.
% Based on Eric Vinande and Dennis Akos
%--------------------------------------------------------------------------
%This program is free software; you can redistribute it and/or
%modify it under the terms of the GNU General Public License
%as published by the Free Software Foundation; either version 2
%of the License, or (at your option) any later version.
%
%This program is distributed in the hope that it will be useful,
%but WITHOUT ANY WARRANTY; without even the implied warranty of
%MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%GNU General Public License for more details.
%
%You should have received a copy of the GNU General Public License
%along with this program; if not, write to the Free Software
%Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301,
%USA.
%--------------------------------------------------------------------------


%% Initialization =========================================================

% C/A code frequency
chipRate        = settings.codeFreqBasis;
% C/A code length
codeLength      = settings.codeLength;
% No. of code periods for coherent integration
cohCodePeriods  = 10;%settings.acquisition.cohCodePeriods;
% Doppler Search Band in Hz
doppSearchBand  = 14*1000;%settings.acqSearchBand*1000;%In Hz
% Original Sampling Frequency
samplingFreq    = settings.samplingFreq;
% Original Sampling period
ts              = 1 / settings.samplingFreq;
% FFT Length should be at least 2*code_length*coh_code_periods
fftLength       = 2^ceil(log2(2*codeLength*cohCodePeriods));
% fftLength = 10*codeLength*cohCodePeriods;
fftLength = samplingFreq/100;
% New Sampling Frequency
reSamplingFreq  = (fftLength/cohCodePeriods)*1000;
% New Sampling period
samplesPerChip  = fftLength / (codeLength * cohCodePeriods);
freqStep        = chipRate / (codeLength * cohCodePeriods);
% No. of frequency bins to search on each side
freqBins        = doppSearchBand/freqStep;
% Samples per each coherent sum
NoOfSamples     = (samplingFreq/chipRate) * codeLength * cohCodePeriods;


%% Remove Carrier and Resample ============================================
% No. of non-coherent summations
nonCohSums      = 10;%settings.acquisition.nonCohSums;


inputfftA = zeros(nonCohSums,fftLength);
inputfftB = zeros(nonCohSums,fftLength);

% Resample parameters
% function rat does not compute an appropriate coefficent when the
% non-coherent period is as large as 300ms.
if rem(cohCodePeriods,2)~=0 && cohCodePeriods >1
    % cohCodePeriods is a odd number
    [resampleNum,resampleDenom] = rat(fftLength/NoOfSamples);
else
    % cohCodePeriods is an even number, except for 1
    commonDiv=gcd(samplingFreq,reSamplingFreq);
    % If the sampling frequency is a prime number
    if commonDiv~=1
        resampleNum=reSamplingFreq/commonDiv;
        resampleDenom=samplingFreq/commonDiv;
    else
        [resampleNum,resampleDenom] = rat(fftLength/NoOfSamples);
    end    
end

% Phase Points
expPhasePoints = exp(i*2*pi*settings.IF*ts*(0:(ceil(2*nonCohSums*NoOfSamples)-1)));
% Remove the carrier and Resample the baseband signal
signalResampled=resample(signal(1:ceil(2*nonCohSums*NoOfSamples)).*expPhasePoints,resampleNum,resampleDenom);
% signalResampled=signal;
for index = 1:nonCohSums

    % Read alternate signal block
    inputA_resamp = signalResampled(2*(index-1)*fftLength+1:2*(index-1)*fftLength+fftLength);
    inputB_resamp = signalResampled((2*index-1)*fftLength+1:2*index*fftLength);

    % Convert the baseband signal to frequency domain
    inputfftA(index,:) = fft(inputA_resamp(1:fftLength));
    inputfftB(index,:) = fft(inputB_resamp(1:fftLength));

end

%--- Initialize acqResults ------------------------------------------------
% Carrier frequencies of detected signals
acqResults.carrFreq     = zeros(1, 158);
% Raw C/A code phases of detected signals
acqResults.rawCodePhase    = zeros(1, 158);
% Refined C/A code phases of detected signals
acqResults.codePhase    = zeros(1, 158);
% Chip rate of C/A code plus doppler
acqResults.codeRate     = zeros(1, 158);
% Correlation peak ratios of the detected signals
acqResults.peakMetric   = zeros(1, 158);

%% Remove Code and Correlate ==============================================
fprintf('(');

%-----For each PRN in the acqSatellite list--------------------------------
for PRN = settings.acqSatelliteList

    %Get the sampled code and compute the FFT
    sampledCode = sampledCAcode(PRN,samplesPerChip*chipRate,fftLength);
    fftCode = fft(sampledCode);
    fftConjCode = conj(fftCode);

    maxCorrAmount = 0;
    freqIndex = 0;

    % For each Doppler search bin
    for binIndex =-freqBins/2:freqBins/2
        
        freqIndex = freqIndex + 1;
        corrA = zeros(1,fftLength);
        corrB = zeros(1,fftLength);

        %Loop through the non-coherent sums
        for index = 1:nonCohSums

            %Circular shift the baseband in frequency domain
            inputfftA_shift = circshift(inputfftA(index,:),[0,binIndex]);
            inputfftB_shift = circshift(inputfftB(index,:),[0,binIndex]);

            %Correlate with the code
            mixed_fftsA = inputfftA_shift .* fftConjCode;
            mixed_fftsB = inputfftB_shift .* fftConjCode;

            %Compute IFFT and sum the results for non-coherent summation
            corrA = corrA + abs(ifft(mixed_fftsA));
            corrB = corrB + abs(ifft(mixed_fftsB));
        end

        %---Find the maximum correlation and store the results---------
        maxCorrA = max(corrA);
        maxCorrB = max(corrB);

        maxCorrCurrentFreq = max(maxCorrA,maxCorrB);

        if ( maxCorrCurrentFreq > maxCorrAmount )
            maxCorrAmount = maxCorrCurrentFreq;
            maxCorrFreqIndex = freqIndex;
            if ( maxCorrB > maxCorrA )
                maxCorr1ms = corrB(1:round(samplesPerChip*codeLength));
            else
                maxCorr1ms = corrA(1:round(samplesPerChip*codeLength));
            end
        end

    end

    acqCarrFreq = (maxCorrFreqIndex - (freqBins/2+1))*freqStep + settings.IF; 

    % Calculate the Peak Height
    [peakHeight,codephaseIndex] = max(maxCorr1ms);
    
    %--- Find 1 chip wide C/A code phase exclude range around the peak ----
    excludeRangeIndex1 = codephaseIndex - ceil(samplesPerChip);
    excludeRangeIndex2 = codephaseIndex + ceil(samplesPerChip);
    
    %--- Correct C/A code phase exclude range if the range includes array
    %boundaries
    if excludeRangeIndex1 < 1  % peak near beginning of code period
        codePhaseRange = excludeRangeIndex2 : ...
            (length(maxCorr1ms) + excludeRangeIndex1);

    elseif excludeRangeIndex2 > length(maxCorr1ms)  % peak near end of code period
        codePhaseRange = (excludeRangeIndex2 - length(maxCorr1ms)) : ...
            excludeRangeIndex1;
    else
        codePhaseRange = [1:excludeRangeIndex1, ...
            excludeRangeIndex2 : length(maxCorr1ms)];
    end

    %--- Find the second highest correlation peak in the same freq. bin ---
    noise_max = max(maxCorr1ms(codePhaseRange));
    acqResults.noise_max(PRN)=noise_max;
    acqResults.maxCorr1ms(PRN,1:length(maxCorr1ms))=maxCorr1ms;
    
    % GPS Satellites
    % Calculate the Peak Metric and the Code Phase
    acqResults.peakMetric(PRN)  = peakHeight / noise_max;
    acqResults.rawCodePhase(PRN)   = round(codephaseIndex/reSamplingFreq*samplingFreq);
    
    
    if (acqResults.peakMetric(PRN)) > settings.acqThreshold
        fprintf('%02d ', PRN);
        
        %Refine the code phase corrected by Xiaofan Li on 9/29/2009
        %Phase Points
        index = refineCodePhase(maxCorr1ms);
        refCodePhase=index/reSamplingFreq*samplingFreq+1;
        acqResults.codePhase(PRN)=round(refCodePhase);
        if (acqResults.codePhase(PRN)==0)  %fix in case it starts at the beginning
            acqResults.codePhase(PRN)=1;
        end

        %Fine Frequency Search only if the coherent integration period is
        %less than 10ms
        if cohCodePeriods < 10
            %signal0DC = signal - mean(signal);
            signal0DC = signal;
            % If the signal is complex, then reform the spectrum back, refer
            % line 108 in postProcessing.m
            if settings.fileType == 2
                signal0DC = conj(signal0DC);
            end
            
            numOfSamples    = round(10*(samplingFreq/chipRate) * codeLength);
            
            %--- Generate 10msec long C/A codes sequence for given PRN --------
            longCaCode = sampledCAcode(PRN,settings.samplingFreq,numOfSamples);
            
            %--- Remove C/A code modulation from the original signal ----------
            
            xCarrierA = ...
                signal0DC(acqResults.codePhase(PRN):(acqResults.codePhase(PRN) + numOfSamples-1)) ...
                .* longCaCode;
            xCarrierB = ...
                signal0DC(acqResults.codePhase(PRN)+numOfSamples:...
                (acqResults.codePhase(PRN) + 2*numOfSamples-1)).* longCaCode;
            
            %--- Find the next highest power of two and increase by 8x --------
            fftNumPts   = 8*(2^(nextpow2(numOfSamples)));
            fftFreqBins = -samplingFreq/2:samplingFreq/fftNumPts:...
                samplingFreq/2-samplingFreq/fftNumPts;
            
            %--- Compute the magnitude of the FFT, find maximum and the
            %associated carrier frequency
            fftxCarrierA=fftshift(abs(fft(xCarrierA, fftNumPts)));
            fftxCarrierB=fftshift(abs(fft(xCarrierB, fftNumPts)));
            
            
            if settings.fileType == 1
                % If it is the real data, then only consider the positive frequency
                % spectrum
                AmpA=abs(fftxCarrierA(fftNumPts/2+1:end));
                AmpB=abs(fftxCarrierB(fftNumPts/2+1:end));
                [vA,maxIndA] = max(AmpA);
                [vB,maxIndB] = max(AmpB);
                maxIndA      = maxIndA+fftNumPts/2;
                maxIndB      = maxIndB+fftNumPts/2;
            else
                % If it is the complex data, consider the whole frequency
                % sepctrum
                AmpA=abs(fftxCarrierA);
                AmpB=abs(fftxCarrierB);
                [vA,maxIndA] = max(AmpA);
                [vB,maxIndB] = max(AmpB);
            end
            
            if vA>vB
                acqResults.carrFreq(PRN)  = fftFreqBins(maxIndA);
            else
                acqResults.carrFreq(PRN)  = fftFreqBins(maxIndB);
            end
            
            % For extreme weak signal the postFFT fine frequency search might not
            % work.
            if abs(acqResults.carrFreq(PRN)-acqCarrFreq) > 750
                acqResults.carrFreq(PRN) = acqCarrFreq;
                %disp('the signal is too weak and postFFT fine frequency search does not work');
            end
        else
            % Do not use postFFT fine frequency search if the coherent
            % integration period is 10ms
            acqResults.carrFreq(PRN) = acqCarrFreq;      
            %disp('Used 10ms coherent integration, no need to use postFFT fine frequency search');
        end
  
        % calculate the doppler of the code based on the doppler of the
        % carrier

        if settings.IF ~= 0
         if settings.carrPhaseInversion == 0
             dopplerCarr=acqResults.carrFreq(PRN)-settings.IF;
        else
             dopplerCarr=settings.IF-acqResults.carrFreq(PRN);
         end
        else
            dopplerCarr=acqResults.carrFreq(PRN);
        end
        dopplerCode=dopplerCarr*settings.codeFreqBasis/1575.42e6;
        acqResults.codeRate(PRN)=settings.codeFreqBasis+dopplerCode;
        
    else
        fprintf('. ');
    end

end


fprintf(')\n');

% save('navSolutions.mat','navSolutions');  