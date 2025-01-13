function [trackResults, channel]= tracking(fid, channel, settings)
% Performs code and carrier tracking for all channels.
%
%[trackResults, channel] = tracking(fid, channel, settings)
%
%   Inputs:
%       fid             - file identifier of the signal record.
%       channel         - PRN, carrier frequencies and code phases of all
%                       satellites to be tracked (prepared by preRum.m from
%                       acquisition results).
%       settings        - receiver settings.
%   Outputs:
%       trackResults    - tracking results (structure array). Contains
%                       in-phase prompt outputs and absolute spreading
%                       code's starting positions, together with other
%                       observation data from the tracking loops. All are
%                       saved every millisecond.

%--------------------------------------------------------------------------
%                           SoftGNSS v3.0
%
% Copyright (C) Dennis M. Akos
% Written by Darius Plausinaitis and Dennis M. Akos
% Based on code by DMAkos Oct-1999
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

%CVS record:
%$Id: tracking.m,v 1.14.2.31 2006/08/14 11:38:22 dpl Exp $

%% Initialize result structure ============================================
% Edited by Sergio Vicenzo - 06 Jan 2025
codePeriods = settings.msToProcess/settings.DPE_cohInt;  % For GPS one C/A code is one ms
% Channel status
trackResults.status         = '-';      % No tracked signal, or lost lock

% The absolute sample in the record of the C/A code start:
trackResults.absoluteSample = zeros(1, codePeriods);

% Freq of the PRN code:
trackResults.codeFreq       = inf(1, codePeriods);

% Frequency of the tracked carrier wave:
trackResults.carrFreq       = inf(1, codePeriods);

% Outputs from the correlators (In-phase):
trackResults.I_P            = zeros(1, codePeriods);
trackResults.I_E            = zeros(1, codePeriods);
trackResults.I_L            = zeros(1, codePeriods);

% Outputs from the correlators (Quadrature-phase):
trackResults.Q_E            = zeros(1, codePeriods);
trackResults.Q_P            = zeros(1, codePeriods);
trackResults.Q_L            = zeros(1, codePeriods);

% Loop discriminators
trackResults.dllDiscr       = inf(1, codePeriods);
trackResults.dllDiscrFilt   = inf(1, codePeriods);
trackResults.pllDiscr       = inf(1, codePeriods);
trackResults.pllDiscrFilt   = inf(1, codePeriods);

% Remain code and carrier phase
trackResults.remCodePhase       = inf(1, codePeriods);
trackResults.remCarrPhase       = inf(1, codePeriods);

%C/No
trackResults.CNo.VSMValue = ...
    zeros(1,floor(codePeriods/settings.CNo.VSMinterval));
trackResults.CNo.VSMIndex = ...
    zeros(1,floor(codePeriods/settings.CNo.VSMinterval));

%--- Copy initial settings for all channels -------------------------------
trackResults = repmat(trackResults, 1, settings.numberOfChannels);

%% Initialize tracking variables ==========================================

% Start waitbar
hwb = waitbar(0,'Tracking...');

%Adjust the size of the waitbar to insert text
CNoPos=get(hwb,'Position');
set(hwb,'Position',[CNoPos(1),CNoPos(2),CNoPos(3),90],'Visible','on');

if (settings.fileType==1)
    dataAdaptCoeff=1;
else
    dataAdaptCoeff=2;
end

%% Start processing channels ==============================================
for channelNr = 1:settings.numberOfChannels

    %% Re-Initialize tracking variables ===================================
    % First run with 1 ms coherent integration to find bit transition
    % Edited by Sergio Vicenzo - 06 Jan 2025
    codePeriods = settings.msToProcess;    % For GPS one C/A code is one ms

    % Only process if PRN is non zero (acquisition was successful)
    if (channel(channelNr).PRN ~= 0)
        % Save additional information - each channel's tracked PRN
        trackResults(channelNr).PRN     = channel(channelNr).PRN;

        % Move the starting point of processing. Can be used to start the
        % signal processing at any point in the data record (e.g. for long
        % records). In addition skip through that data file to start at the
        % appropriate sample (corresponding to code phase). Assumes sample
        % type is schar (or 1 byte per sample)

        % Edited by Sergio Vicenzo to also handle "int16" files
        % 12 Feb 2024
        if strcmp(settings.dataType,'int16')  
            fseek(fid, ...
                dataAdaptCoeff*(settings.skipNumberOfBytes + ...
                (channel(channelNr).codePhase-1)*2), ...
                'bof');
        else
            fseek(fid, ...
                dataAdaptCoeff*(settings.skipNumberOfBytes + ...
                channel(channelNr).codePhase-1), ...
                'bof');
        end

        % Get a vector with the C/A code sampled 1x/chip
        caCode = generateCAcode(channel(channelNr).PRN);
        % Then make it possible to do early and late versions
        caCode = [caCode(1023) caCode caCode(1)];

        %--- Perform various initializations ------------------------------

        % define initial code frequency basis of NCO
        codeFreq      = settings.codeFreqBasis;
        % define residual code phase (in chips)
        remCodePhase  = 0.0;
        % define carrier frequency which is used over whole tracking period
        carrFreq      = channel(channelNr).acquiredFreq;
        carrFreqBasis = channel(channelNr).acquiredFreq;
        % define residual carrier phase
        remCarrPhase  = 0.0;

        %code tracking loop parameters
        oldCodeNco   = 0.0;
        oldCodeError = 0.0;

        %carrier/Costas loop parameters
        oldCarrNco   = 0.0;
        oldCarrError = 0.0;

        settings.dllNoiseBandwidth       = 1.5;
        settings.pllNoiseBandwidth       = 20;

        %--- DLL variables --------------------------------------------------------
        % Define early-late offset (in chips)
        earlyLateSpc = settings.dllCorrelatorSpacing;
        % Calculate filter coefficient values
        [tau1code, tau2code] = calcLoopCoef(settings.dllNoiseBandwidth, ...
            settings.dllDampingRatio, ...
            1.0);
        % Summation interval
        PDIcode = 0.001;
        
        %--- PLL variables --------------------------------------------------------
        % Calculate filter coefficient values
        [tau1carr, tau2carr] = calcLoopCoef(settings.pllNoiseBandwidth, ...
                                            settings.pllDampingRatio, ...
                                            0.25);
        % Summation interval
        PDIcarr = 0.001;

        %=== Process the number of specified code periods =================
        for loopCnt =  1:codePeriods

            %% GUI update -------------------------------------------------------------
            % The GUI is updated every 50ms. This way Matlab GUI is still
            % responsive enough. At the same time Matlab is not occupied
            % all the time with GUI task.
            % Commented by Sergio Vicenzo - 06 Jan 2025
%             if (rem(loopCnt, 50) == 0)
% 
%                 Ln=sprintf('\n');
%                 trackingStatus=['Tracking: Ch ', int2str(channelNr), ...
%                     ' of ', int2str(settings.numberOfChannels),Ln ...
%                     'PRN: ', int2str(channel(channelNr).PRN),Ln ...
%                     'Completed ',int2str(loopCnt), ...
%                     ' of ', int2str(codePeriods), ' msec',Ln...
%                     'C/No: ',CNo,' (dB-Hz)'];
% 
%                 try
%                     waitbar(loopCnt/codePeriods, ...
%                         hwb, ...
%                         trackingStatus);
%                 catch
%                     % The progress bar was closed. It is used as a signal
%                     % to stop, "cancel" processing. Exit.
%                     disp('Progress bar closed, exiting...');
%                     return
%                 end
%             end

            %% Read next block of data ------------------------------------------------
            % Record sample number (based on 8bit samples)
            % Edited by Sergio Vicenzo to also handle "int16" files
            % 12 Feb 2024
            if strcmp(settings.dataType,'int16') 
                trackResults(channelNr).absoluteSample(loopCnt) = ...
                    (ftell(fid))/dataAdaptCoeff/2;
            else
                trackResults(channelNr).absoluteSample(loopCnt) = ...
                    (ftell(fid))/dataAdaptCoeff;
            end

			% Find the size of a "block" or code period in whole samples
            % Update the phasestep based on code freq (variable) and
            % sampling frequency (fixed)
            codePhaseStep = codeFreq / settings.samplingFreq;
            
            % Find the size of a "block" or code period in whole samples
            blksize = ceil((settings.codeLength-remCodePhase) / codePhaseStep);

            % Read in the appropriate number of samples to process this
            % interation
            [rawSignal, samplesRead] = fread(fid, ...
                dataAdaptCoeff*blksize, settings.dataType);

            rawSignal = rawSignal';

            if (dataAdaptCoeff==2)
                rawSignal1=rawSignal(1:2:end);
                rawSignal2=rawSignal(2:2:end);
                rawSignal = rawSignal1 + 1i .* rawSignal2;  % transpose vector
            end


            % If did not read in enough samples, then could be out of
            % data - better exit
            if (samplesRead ~= dataAdaptCoeff*blksize)
                disp('Not able to read the specified number of samples  for tracking, exiting!')
                fclose(fid);
                return
            end

            %% Set up all the code phase tracking information -------------------------
            % Save remCodePhase for current correlation
            trackResults(channelNr).remCodePhase(loopCnt) = remCodePhase;
            
            % Define index into early code vector
            tcode       = (remCodePhase-earlyLateSpc) : ...
                codePhaseStep : ...
                ((blksize-1)*codePhaseStep+remCodePhase-earlyLateSpc);
            tcode2      = ceil(tcode) + 1;
            earlyCode   = caCode(tcode2);

            % Define index into late code vector
            tcode       = (remCodePhase+earlyLateSpc) : ...
                codePhaseStep : ...
                ((blksize-1)*codePhaseStep+remCodePhase+earlyLateSpc);
            tcode2      = ceil(tcode) + 1;
            lateCode    = caCode(tcode2);

            % Define index into prompt code vector
            tcode       = remCodePhase : ...
                codePhaseStep : ...
                ((blksize-1)*codePhaseStep+remCodePhase);
            tcode2      = ceil(tcode) + 1;
            promptCode  = caCode(tcode2);

            remCodePhase = (tcode(blksize) + codePhaseStep) - settings.codeLength;

            %% Generate the carrier frequency to mix the signal to baseband -----------
            % Save remCarrPhase for current correlation
            trackResults(channelNr).remCarrPhase(loopCnt) = remCarrPhase;

            % Get the argument to sin/cos functions
            time    = (0:blksize) ./ settings.samplingFreq;
            trigarg = ((carrFreq * 2.0 * pi) .* time) + remCarrPhase;
            remCarrPhase = rem(trigarg(blksize+1), (2 * pi));

            % Finally compute the signal to mix the collected data to
            % bandband
            carrsig = exp(1i .* trigarg(1:blksize));

            %% Generate the six standard accumulated values ---------------------------
            % First mix to baseband
            qBasebandSignal = real(carrsig .* rawSignal);
            iBasebandSignal = imag(carrsig .* rawSignal);

            % Now get early, late, and prompt values for each
            I_E = sum(earlyCode  .* iBasebandSignal);
            Q_E = sum(earlyCode  .* qBasebandSignal);
            I_P = sum(promptCode .* iBasebandSignal);
            Q_P = sum(promptCode .* qBasebandSignal);
            I_L = sum(lateCode   .* iBasebandSignal);
            Q_L = sum(lateCode   .* qBasebandSignal);

            %% Find PLL error and update carrier NCO ----------------------------------

            % Implement carrier loop discriminator (phase detector)
            carrError = atan(Q_P / I_P) / (2.0 * pi);

            % Implement carrier loop filter and generate NCO command
            carrNco = oldCarrNco + (tau2carr/tau1carr) * ...
                (carrError - oldCarrError) + carrError * (PDIcarr/tau1carr);
            oldCarrNco   = carrNco;
            oldCarrError = carrError;

            % Save carrier frequency for current correlation
            trackResults(channelNr).carrFreq(loopCnt) = carrFreq;

            % Modify carrier freq based on NCO command
            carrFreq = carrFreqBasis + carrNco;

            

            %% Find DLL error and update code NCO -------------------------------------
            codeError = (sqrt(I_E * I_E + Q_E * Q_E) - sqrt(I_L * I_L + Q_L * Q_L)) / ...
                (sqrt(I_E * I_E + Q_E * Q_E) + sqrt(I_L * I_L + Q_L * Q_L));

            % Implement code loop filter and generate NCO command
            codeNco = oldCodeNco + (tau2code/tau1code) * ...
                (codeError - oldCodeError) + codeError * (PDIcode/tau1code);
            oldCodeNco   = codeNco;
            oldCodeError = codeError;
            
            % Save code frequency for current correlation
            trackResults(channelNr).codeFreq(loopCnt) = codeFreq;

            % Modify code freq based on NCO command
            codeFreq = settings.codeFreqBasis - codeNco;

            %% Record various measures to show in postprocessing ----------------------
            trackResults(channelNr).dllDiscr(loopCnt)       = codeError;
            trackResults(channelNr).dllDiscrFilt(loopCnt)   = codeNco;
            trackResults(channelNr).pllDiscr(loopCnt)       = carrError;
            trackResults(channelNr).pllDiscrFilt(loopCnt)   = carrNco;

            trackResults(channelNr).I_E(loopCnt) = I_E;
            trackResults(channelNr).I_P(loopCnt) = I_P;
            trackResults(channelNr).I_L(loopCnt) = I_L;
            trackResults(channelNr).Q_E(loopCnt) = Q_E;
            trackResults(channelNr).Q_P(loopCnt) = Q_P;
            trackResults(channelNr).Q_L(loopCnt) = Q_L;

            % Commented by Sergio Vicenzo - 07 Jan 2025
%             if (rem(loopCnt,settings.CNo.VSMinterval)==0)
%                 vsmCnt=vsmCnt+1;
%                 CNoValue=CNoVSM(trackResults(channelNr).I_P(loopCnt-settings.CNo.VSMinterval+1:loopCnt),...
%                     trackResults(channelNr).Q_P(loopCnt-settings.CNo.VSMinterval+1:loopCnt),settings.CNo.accTime);
%                 trackResults(channelNr).CNo.VSMValue(vsmCnt)=CNoValue;
%                 trackResults(channelNr).CNo.VSMIndex(vsmCnt)=loopCnt;
%                 CNo=int2str(CNoValue);
%             end

            % Check for bit transition --------------------------------------------
            % Locate bit transition after 200 ms of tracking with PDI = 1
            % ms
            % 200 ms was chosen arbitrarily based on observation that the
            % tracking loop has stabilised
            % Added by Sergio Vicenzo - 06 Jan 2025
            %
            if loopCnt > 200 && ...
                    sign(I_P)~=sign(trackResults(channelNr).I_P(loopCnt-1))

                newStart=loopCnt;
                break;

            end

        end % for loopCnt

        %% Restart Tracking with desired coherent integration time ================
        % Added by Sergio Vicenzo - 06 Jan 2025
        % Re-Initialize tracking variables ----------------------------------------

        settings.dllNoiseBandwidth       = 4;
        settings.pllNoiseBandwidth       = 5;

        % Get a vector with the C/A code sampled 1x/chip
        caCode = generateCAcode(channel(channelNr).PRN);
        % Then make it possible to do early and late versions
        caCode = [caCode(1023) repmat(caCode,1,settings.DPE_cohInt) caCode(1)];

        ori_codeperiods=settings.msToProcess;
        codePeriods=round(codePeriods/settings.DPE_cohInt);

        %--- DLL variables --------------------------------------------------------             
        % Summation interval
        PDIcode = settings.DPE_cohInt/1000;

        % Calculate filter coefficient values
        [tau1code, tau2code] = calcLoopCoef(settings.dllNoiseBandwidth, ...
            settings.dllDampingRatio, ...
            1.0);
        
                
        %--- PLL variables --------------------------------------------------------
        % Summation interval
        PDIcarr = settings.DPE_cohInt/1000;

        % Calculate filter coefficient values
        [tau1carr, tau2carr] = calcLoopCoef(settings.pllNoiseBandwidth, ...
                                            settings.pllDampingRatio, ...
                                            0.25);

        %--- Perform various initializations --------------------------------------

        % define initial code frequency basis of NCO
        codeFreq      = trackResults(channelNr).codeFreq(newStart-1);
        % define residual code phase (in chips)
        remCodePhase  = trackResults(channelNr).remCodePhase(newStart);
        % define carrier frequency which is used over whole tracking period
        carrFreq      = trackResults(channelNr).carrFreq(newStart-1);
        carrFreqBasis = channel(channelNr).acquiredFreq;
        % define residual carrier phase
        remCarrPhase  = trackResults(channelNr).remCarrPhase(newStart);

        %code tracking loop parameters
        oldCodeNco   = trackResults(channelNr).dllDiscrFilt(newStart-1);
        oldCodeError = trackResults(channelNr).dllDiscr(newStart-1);

        %carrier/Costas loop parameters
        oldCarrNco   = trackResults(channelNr).pllDiscrFilt(newStart-1);
        oldCarrError = trackResults(channelNr).pllDiscr(newStart-1);

        %C/No computation
        vsmCnt  = 0;CNo = 0;  
        
        %=== Process the number of specified code periods =================
        for loopCnt =  1:codePeriods

            %% GUI update -------------------------------------------------------------
            % The GUI is updated every 50ms. This way Matlab GUI is still
            % responsive enough. At the same time Matlab is not occupied
            % all the time with GUI task.
            if (rem(loopCnt, 50) == 0)

                Ln=sprintf('\n');
                trackingStatus=['Tracking: Ch ', int2str(channelNr), ...
                    ' of ', int2str(settings.numberOfChannels),Ln ...
                    'PRN: ', int2str(channel(channelNr).PRN),Ln ...
                    'Completed ',int2str(loopCnt* settings.DPE_cohInt), ...
                    ' of ', int2str(ori_codeperiods), ' msec',Ln...
                    'C/No: ',CNo,' (dB-Hz)'];

                try
                    waitbar(loopCnt/codePeriods, ...
                        hwb, ...
                        trackingStatus);
                catch
                    % The progress bar was closed. It is used as a signal
                    % to stop, "cancel" processing. Exit.
                    disp('Progress bar closed, exiting...');
                    return
                end
            end

            % Move the starting point of processing. Can be used to start the
            % signal processing at any point in the data record (e.g. for long
            % records). In addition skip through that data file to start at the
            % appropriate sample (corresponding to code phase). Assumes sample
            % type is schar (or 1 byte per sample)
    
            % Edited by Sergio Vicenzo to also handle "int16" files
            % 06 Jan 2025
            if loopCnt==1
                if strcmp(settings.dataType,'int16') 
                    fseek(fid, ...
                        dataAdaptCoeff*(...
                        (trackResults(channelNr).absoluteSample(newStart))*2), ...
                        'bof');
                else
                    fseek(fid, ...
                        dataAdaptCoeff*(trackResults(channelNr).absoluteSample(newStart)), ...
                        'bof');
                end
            end

            %% Read next block of data ------------------------------------------------
            % Record sample number (based on 8bit samples)
            % Edited by Sergio Vicenzo to also handle "int16" files
            % 12 Feb 2024
            if strcmp(settings.dataType,'int16') 
                trackResults(channelNr).absoluteSample(loopCnt) = ...
                    (ftell(fid))/dataAdaptCoeff/2;
            else
                trackResults(channelNr).absoluteSample(loopCnt) = ...
                    (ftell(fid))/dataAdaptCoeff;
            end

			% Find the size of a "block" or code period in whole samples
            % Update the phasestep based on code freq (variable) and
            % sampling frequency (fixed)

            codePhaseStep = codeFreq / settings.samplingFreq;

            
            % Find the size of a "block" or code period in whole samples

            blksize = ceil(((settings.codeLength*settings.DPE_cohInt)-remCodePhase) ...
            / codePhaseStep);



            % Read in the appropriate number of samples to process this
            % interation
            [rawSignal, samplesRead] = fread(fid, ...
                dataAdaptCoeff*blksize, settings.dataType);

            rawSignal = rawSignal';

            if (dataAdaptCoeff==2)
                rawSignal1=rawSignal(1:2:end);
                rawSignal2=rawSignal(2:2:end);
                rawSignal = rawSignal1 + 1i .* rawSignal2;  % transpose vector
            end


            % If did not read in enough samples, then could be out of
            % data - better exit
            if (samplesRead ~= dataAdaptCoeff*blksize)
                disp('Not able to read the specified number of samples  for tracking, exiting!')
                fclose(fid);
                return
            end

            %% Set up all the code phase tracking information -------------------------
            % Save remCodePhase for current correlation

            trackResults(channelNr).remCodePhase(loopCnt) = remCodePhase;


            
            % Define index into early code vector
            tcode       = (remCodePhase-earlyLateSpc) : ...
                codePhaseStep : ...
                ((blksize-1)*codePhaseStep+remCodePhase-earlyLateSpc);
            tcode2      = ceil(tcode) + 1;
            earlyCode   = caCode(tcode2);

            % Define index into late code vector
            tcode       = (remCodePhase+earlyLateSpc) : ...
                codePhaseStep : ...
                ((blksize-1)*codePhaseStep+remCodePhase+earlyLateSpc);
            tcode2      = ceil(tcode) + 1;
            lateCode    = caCode(tcode2);

            % Define index into prompt code vector
            tcode       = remCodePhase : ...
                codePhaseStep : ...
                ((blksize-1)*codePhaseStep+remCodePhase);
            tcode2      = ceil(tcode) + 1;
            promptCode  = caCode(tcode2);

            remCodePhase = (tcode(blksize) + codePhaseStep)...
                - settings.codeLength * (settings.DPE_cohInt);

            %% Generate the carrier frequency to mix the signal to baseband -----------
            % Save remCarrPhase for current correlation

            trackResults(channelNr).remCarrPhase(loopCnt) = remCarrPhase;


            % Get the argument to sin/cos functions
            time    = (0:blksize) ./ settings.samplingFreq;
            trigarg = ((carrFreq * 2.0 * pi) .* time) + remCarrPhase;
            remCarrPhase = rem(trigarg(blksize+1), (2 * pi));

            % Finally compute the signal to mix the collected data to
            % bandband
            carrsig = exp(1i .* trigarg(1:blksize));

            %% Generate the six standard accumulated values ---------------------------
            % First mix to baseband
            qBasebandSignal = real(carrsig .* rawSignal);
            iBasebandSignal = imag(carrsig .* rawSignal);

            % Now get early, late, and prompt values for each
            I_E = sum(earlyCode  .* iBasebandSignal);
            Q_E = sum(earlyCode  .* qBasebandSignal);
            I_P = sum(promptCode .* iBasebandSignal);
            Q_P = sum(promptCode .* qBasebandSignal);
            I_L = sum(lateCode   .* iBasebandSignal);
            Q_L = sum(lateCode   .* qBasebandSignal);

            %% Find PLL error and update carrier NCO ----------------------------------

            % Implement carrier loop discriminator (phase detector)
            carrError = atan(Q_P / I_P) / (2.0 * pi);

            % Implement carrier loop filter and generate NCO command

                
            carrNco = oldCarrNco + (tau2carr/tau1carr) * ...
                (carrError - oldCarrError) + carrError * (PDIcarr/tau1carr);
            oldCarrNco   = carrNco;
            oldCarrError = carrError;

            % Save carrier frequency for current correlation
            trackResults(channelNr).carrFreq(loopCnt) = carrFreq;

            % Modify carrier freq based on NCO command
            carrFreq = carrFreqBasis + carrNco;

            

            %% Find DLL error and update code NCO -------------------------------------
            codeError = (sqrt(I_E * I_E + Q_E * Q_E) - sqrt(I_L * I_L + Q_L * Q_L)) / ...
                (sqrt(I_E * I_E + Q_E * Q_E) + sqrt(I_L * I_L + Q_L * Q_L));

            % Implement code loop filter and generate NCO command

            codeNco = oldCodeNco + (tau2code/tau1code) * ...
            (codeError - oldCodeError) + codeError * (PDIcode/tau1code);
            oldCodeNco   = codeNco;
            oldCodeError = codeError;
            
            % Save code frequency for current correlation
            trackResults(channelNr).codeFreq(loopCnt) = codeFreq;

            % Modify code freq based on NCO command
            codeFreq = settings.codeFreqBasis - codeNco;

            %% Record various measures to show in postprocessing ----------------------
            trackResults(channelNr).dllDiscr(loopCnt)       = codeError;
            trackResults(channelNr).dllDiscrFilt(loopCnt)   = codeNco;
            trackResults(channelNr).pllDiscr(loopCnt)       = carrError;
            trackResults(channelNr).pllDiscrFilt(loopCnt)   = carrNco;

            trackResults(channelNr).I_E(loopCnt) = I_E;
            trackResults(channelNr).I_P(loopCnt) = I_P;
            trackResults(channelNr).I_L(loopCnt) = I_L;
            trackResults(channelNr).Q_E(loopCnt) = Q_E;
            trackResults(channelNr).Q_P(loopCnt) = Q_P;
            trackResults(channelNr).Q_L(loopCnt) = Q_L;

            if (rem(loopCnt,settings.CNo.VSMinterval)==0)
                vsmCnt=vsmCnt+1;
                CNoValue=CNoVSM(trackResults(channelNr).I_P(loopCnt-settings.CNo.VSMinterval+1:loopCnt),...
                    trackResults(channelNr).Q_P(loopCnt-settings.CNo.VSMinterval+1:loopCnt),...
                    settings.CNo.accTime*settings.DPE_cohInt);
                trackResults(channelNr).CNo.VSMValue(vsmCnt)=CNoValue;
                trackResults(channelNr).CNo.VSMIndex(vsmCnt)=loopCnt;
                CNo=int2str(CNoValue);
            end

          

        end % for loopCnt

        % If we got so far, this means that the tracking was successful
        % Now we only copy status, but it can be update by a lock detector
        % if implemented
        trackResults(channelNr).status  = channel(channelNr).status;

    end % if a PRN is assigned
end % for channelNr

% Close the waitbar
close(hwb)
