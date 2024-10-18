function [navSolutions]=DPE_module...
    (currMeasNr,navSolutions,activeChnList,trackResults,currMeasSample,...
    satPositions,transmitTime,localTime,codePhase,settings,satElev,fid,...
    LS_clkBias,satClkCorr)
%
% DPE_module is a Direct Position Estimation (DPE) plug-in module that can 
% be integrated into existing two-step positioning MATLAB SDRs.
%
% The correlations per every code phase would first be computed and later
% interpolated to obtain correlations for every candidate positions
% (a correlogram over the candidate positions)
%
% The plug-in module is set by default to use 1 ms for both coherent and 
% non-coherent integration
%
% [navSolutions]=DPE_module...
%   (currMeasNr,navSolutions,activeChnList,trackResults,currMeasSample,...
%   satPositions, transmitTime,localTime,codePhase,settings,satElev,fid,...
%   LS_clkBias,satClkCorr)

%   Inputs:
%       currMeasNr      - Current positioning epoch
%       navSolutions    - Contains Least Squares positioning estimates 
%                         for DPE's search space initialization
%       activeChnList   - List of channels to be processed
%       trackResults    - Output from the tracking function
%       currMeasSample  - Current measurement sample location
%                         (measurement time)
%       satPositions    - Satellites positions (in ECEF system: [X; Y; Z;] -
%                         one column per satellite) from Least Squares
%       transmitTime    - Transmitting time of channels to be processed 
%                         corresponding to measurement time 
%       localTime       - Local time(in GPST) at measurement time,
%                         uncorrected by Least Square's clock bias
%       codePhase       - Code phase from start of a PRN code to current
%                         measurement sample location of all satellites
%       settings        - Receiver settings 
%       satElev         - Satellites elevation angles (degrees) 
%                         from Least Squares
%       fid             - Integer file identifier for raw IF data
%       LS_clkBias      - Receiver clock bias estimated by Least Squares
%       satClkCorr      - Satellite clock bias

%
%   Outputs:
%       DPE_latitude    - Estimated receiver latitude from DPE
%       DPE_longitude   - Estimated receiver longitude from DPE
%       DPE_height      - Estimated receiver height from DPE
%       DPE_clkBias     - Estimated receiver clock bias from DPE

%--------------------------------------------------------------------------
%------------------------- DPE_module v1.0 --------------------------------
%   Copyright (C) 2024 Sergio Vicenzo
%   Written by Sergio Vicenzo
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
%   This program is free software; you can redistribute it and/or modify
%   it under the terms of the GNU General Public License as published by
%   the Free Software Foundation; either version 2 of the License, or
%   (at your option) any later version.
% 
%   This program is distributed in the hope that it will be useful,
%   but WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%   GNU General Public License for more details.
% 
%   You should have received a copy of the GNU General Public License along
%   with this program; if not, write to the Free Software Foundation, Inc.,
%   51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.

% Last Updated: 7 June 2024

addpath include common

% === Initialize parameters ===============================================
dtr                     = pi/180;
m2lat                   = 1/110734;
m2lon                   = 1/103043;
trop                    = zeros(length(activeChnList),1);
chip_spacings           = [-(flip(settings.chipspacing_dpe_precalc:...
                            settings.chipspacing_dpe_precalc:1)),0,...
                            settings.chipspacing_dpe_precalc:...
                            settings.chipspacing_dpe_precalc:1];
precalc_correlations    = ...
                          zeros(length(chip_spacings)...
                          *length(activeChnList),3);

% === Make folder to store correlograms ===================================
if ~exist([settings.outfile_root '\Correlogram\'], 'dir') && ...
        settings.DPE_plotCorrelogram == 1
    mkdir([settings.outfile_root '\Correlogram\']);
end

% === Generate Latitude and Longitude search space ========================
% Lat-long search space is of plus-minus "settings.DPE_latlong_span" meters
% Centered at Least Squares' lat-long estimate
[~,~,~,mesh_ll] = ...
   meshgrid_LIST(navSolutions.latitude(currMeasNr)+(m2lat).*...
   (-settings.DPE_latlong_span:settings.candPVT_spacing:settings.DPE_latlong_span)...
    ,navSolutions.longitude(currMeasNr)+(m2lon).*...
    (-settings.DPE_latlong_span:settings.candPVT_spacing:settings.DPE_latlong_span));

% === Generate Height search space ========================================
% Height search space is of plus-minus "settings.DPE_height_span" meters
% Centered at Least Squares' height estimate
alt_search = round(navSolutions.height(currMeasNr),0)-...
    settings.DPE_height_span:settings.candPVT_spacing:...
    round(navSolutions.height(currMeasNr),0)+settings.DPE_height_span;

% === Generate Clock Bias search space ====================================
% Use a wider search space for first epoch to ensure clock bias
% convergence 
if currMeasNr == 1
    candia_xyz_LS = [navSolutions.X(currMeasNr),...
        navSolutions.Y(currMeasNr),navSolutions.Z(currMeasNr)];
    clock_bias=zeros(length(activeChnList),1);
    simulated_pseudorange = zeros(length(activeChnList),1);
    for i = 1:length(activeChnList)
        if (settings.useTropCorr == 1)
         % === Estimate tropospheric error for each satellite =============
         trop(i) = tropo(sin(satElev(activeChnList(i)) * dtr), ...
             0.0, 1013.0, 293.0, 50.0, 0.0, 0.0, 0.0);
        end
    
        simulated_pseudorange(i)=sqrt(sum((satPositions(:,i)'...
        - candia_xyz_LS).^2));
    
        clock_bias(i) = ...
            (navSolutions.rawP(activeChnList(i),currMeasNr)+ ...
               satClkCorr(i) * settings.c - simulated_pseudorange(i));
    end

    clk_s_area=...
        round(min(clock_bias))-100:1:round(max(clock_bias))+100;

else

% Since Least Squares' clock bias is already stable and relatively accurate
% after the first epoch, use narrow search space in consecutive positioning 
% epochs
    clk_s_area=round((LS_clkBias))-settings.DPE_clkBias_span:1:...
        round((LS_clkBias))+settings.DPE_clkBias_span;

    if (settings.useTropCorr == 1)
        for i = 1:length(activeChnList)
         % === Estimate tropospheric error for each satellite =============
         trop(i) = tropo(sin(satElev(activeChnList(i)) * dtr), ...
             0.0, 1013.0, 293.0, 50.0, 0.0, 0.0, 0.0);
        end
    end
end

% If clock bias search space is too large
if length(clk_s_area)>100000 % remember to uncomment this
     clk_s_area=round((LS_clkBias))-settings.DPE_clkBias_span:1:...
        round((LS_clkBias))+settings.DPE_clkBias_span;
end

% === Pre-calculate the correlation values ================================
% A pure, traditional DPE implementation would require iterative
% correlations for every candidate position, which is highly
% computationally extensive. This section pre-calculates the correlations
% in order of "settings.chipspacing_dpe_precalc" chip spacing
% to save computational time

for j=1:length(activeChnList)
       
        % === Obtain closest tracking measurements to the current =========
        % === measurement sample ==========================================

        [~,closestIndex] = min(abs(currMeasSample-...
            trackResults(activeChnList(j)).absoluteSample));
        if currMeasSample < ...
              trackResults(activeChnList(j)).absoluteSample(closestIndex)

          closestIndex=closestIndex-1;
        end

        % === Move IF data to the current measurement sample ==============
        % fread would move the measurement sample automatically for every 
        % iteration of j. This code prevents the current sample to move.
        if closestIndex > (settings.DPE_nonCohInt-1)
            nonCohInt=closestIndex-settings.DPE_nonCohInt+1:closestIndex;
        else
            nonCohInt=1:closestIndex;
        end

        % Generate indexes to store the pre-calculated correlations in an
        % array
        count2=3;
        count=1+(j-1)*length(chip_spacings);

     % === Pre-calculate the correlations =================================
     for closestIndex=nonCohInt
        % Loop repeats depending on the non-coherent integration time
        % chosen

        % === Record sample number (based on 8-bit and 16-bit samples) ====
        if strcmp(settings.dataType,'int16') 
            fseek(fid, ...
                settings.fileType*...
                ((trackResults(activeChnList(j)).absoluteSample(closestIndex)-1)*2),'bof');
        else
            fseek(fid, ...
                settings.fileType*...
                (trackResults(activeChnList(j)).absoluteSample(closestIndex)-1),'bof');
        end

        % === Update the phasestep based on code freq (variable) and ======
        % === sampling frequency (fixed) ==================================
        codePhaseStep = ...
            trackResults(activeChnList(j)).codeFreq(closestIndex)...
            / settings.samplingFreq;
               
        % === Find the size of a "block" or code period in whole samples ==
        blksize = ceil((settings.codeLength-...
            trackResults(activeChnList(j)).remCodePhase(closestIndex)) ...
            / codePhaseStep);

        % === Read in the appropriate number of samples to process =======
        [rawSignal, ~] = fread(fid, ...
           settings.fileType*blksize, settings.dataType);

        rawSignal = rawSignal';

        if (settings.fileType==2)
           rawSignal1=rawSignal(1:2:end);
           rawSignal2=rawSignal(2:2:end);
           rawSignal = rawSignal1 + 1i .* rawSignal2;%transpose vector
        end
            
        % === Get the argument to sin/cos functions =======================
        time    = (0:blksize) ./ settings.samplingFreq;
        remCarrPhase = ...
            trackResults(activeChnList(j)).remCarrPhase(closestIndex);
        trigarg = ...
            ((trackResults(activeChnList(j)).carrFreq(closestIndex)...
            * 2.0 * pi) .* time) + remCarrPhase;

        % === Compute the carrier signal ==================================
        carrsig = exp(1i .* trigarg(1:blksize));

        % === Mix the carrier replica to baseband =========================
        qBasebandSignal = real(carrsig .* rawSignal);
        iBasebandSignal = imag(carrsig .* rawSignal);

        % === Get the C/A code ============================================
        caCode = generateCAcode(trackResults(activeChnList(j)).PRN);

        for spacing=round(chip_spacings,2)

        % === Define index into the code vector ===========================
        delay_index = ...
            trackResults(activeChnList(j)).remCodePhase(closestIndex)-spacing: ...
            codePhaseStep : ...
            ((blksize-1)*codePhaseStep + ...
            trackResults(activeChnList(j)).remCodePhase(closestIndex)-spacing);       

            caCode1 = [caCode(end) caCode caCode(1)];
            tcodee = ceil(delay_index)+1;
            tcodee(tcodee==0)=tcodee(tcodee==0)+1;

            s = caCode1(tcodee);

            I = sum(s  .* iBasebandSignal);
            Q = sum(s  .* qBasebandSignal);

        % === Define corresponding code phase of the computed correlation =
        if codePhase(activeChnList(j))+spacing < 0
            codePhase5 = ...
                codePhase(activeChnList(j))+spacing +settings.codeLength;
        elseif codePhase(activeChnList(j))+spacing >settings.codeLength
            codePhase5 = ...
                abs(settings.codeLength-(codePhase(activeChnList(j))+spacing));
        else
            codePhase5 = codePhase(activeChnList(j))+spacing;
        end
            
        % === Store the correlations and its corresponding code phase =====
          precalc_correlations(count,1) = ...
              trackResults(activeChnList(j)).PRN;
          precalc_correlations(count,2) = codePhase5;
          precalc_correlations(count,count2) = sqrt(I.^2 + Q.^2);

          count=count+1;
        end
          count=1+(j-1)*length(chip_spacings);
          count2=count2+1;
     end 
     % End for pre-calculating the corr values for a single satellite

end % End for pre-calculating the corr values for all satellites

% === Non-coherent integration of each satellite's autocorrelation function
precalc_correlations(:,3) = sum(precalc_correlations(1:end,3:end),2);


% X_var = optimizableVariable('X',[min(candia_xyz_clkBias(:,1)) max(candia_xyz_clkBias(:,1))]);
% Y_var = optimizableVariable('Y',[min(candia_xyz_clkBias(:,2)) max(candia_xyz_clkBias(:,2))]);
% Z_var = optimizableVariable('Z',[min(candia_xyz_clkBias(:,3)) max(candia_xyz_clkBias(:,3))]);

X_var = optimizableVariable('X',[min(mesh_ll(:,1)) max(mesh_ll(:,1))]);
Y_var = optimizableVariable('Y',[min(mesh_ll(:,2)) max(mesh_ll(:,2))]);
Z_var = optimizableVariable('Z',[min(alt_search) max(alt_search)]);
clkBias_var = optimizableVariable('clkBias',[min(clk_s_area) max(clk_s_area)]);

fun = @(x) -computeCorrelation([x.X x.Y x.Z],...
                satPositions,trop,satClkCorr,x.clkBias,localTime,...
                transmitTime,codePhase,precalc_correlations,...
                activeChnList,trackResults,settings);

results = bayesopt(fun,[X_var,Y_var,Z_var,clkBias_var],...
    'MaxObjectiveEvaluations',80,...
    'NumSeedPoints',20,'MaxTime',100,...
    'Verbose',0,'PlotFcn',[],'AcquisitionFunctionName', 'expected-improvement');%,...
%     'OutputFcn',@terminationFunc);%{@plotObjectiveModel,@plotMinObjective});%,'UseParallel',true);%,'InitialX', ...
%     table(init_candia_llh_clkBias(1), init_candia_llh_clkBias(2), init_candia_llh_clkBias(3), init_candia_llh_clkBias(4)),...
%     'InitialObjective',init_corr_val);


% === Record DPE estimates ================================================
navSolutions.DPE_estimate(currMeasNr, 1:5) =  ...
    [results.XAtMinObjective.X...
    results.XAtMinObjective.Y...
    results.XAtMinObjective.Z...
    -results.MinObjective,...
    results.XAtMinObjective.clkBias];

navSolutions.DPE_latitude(currMeasNr) =  ...
    results.XAtMinObjective.X;

navSolutions.DPE_longitude(currMeasNr) =  ...
    results.XAtMinObjective.Y;

navSolutions.DPE_height(currMeasNr) =  ...
    results.XAtMinObjective.Z;

navSolutions.DPE_clkBias(currMeasNr) = ...
    results.XAtMinObjective.clkBias;


close all
end


