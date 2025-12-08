function [navSolutions]=CorrDPE_ARS...
    (currMeasNr,navSolutions,activeChnList,trackResults,currMeasSample,...
    satPositions,transmitTime,localTime,settings,satElev,...
    LS_clkBias,satClkCorr)

% CorrDPE_ARS is a Corr-DPE plug-in module that can be integrated into 
% existing two-step positioning MATLAB SDRs.
%
% Accelerated Random Search (ARS) is employed for Corr-DPE to obtain the 
% PVT estimates.
%
% [navSolutions]= CorrDPE_ARS...
%   (currMeasNr,navSolutions,activeChnList,trackResults,currMeasSample,...
%   satPositions, transmitTime,localTime,codePhase,settings,satElev,...
%   LS_clkBias,satClkCorr)

%   Inputs:
%       currMeasNr      - Current positioning epoch
%       navSolutions    - Contains Least Squares positioning estimates 
%                         for DPE's search space initialization
%       activeChnList   - List of channels to be processed
%       trackResults    - Output from the tracking function, contains the
%                         multicorrelator values to be used for positioning
%       currMeasSample  - Current measurement sample location
%                         (measurement time)
%       satPositions    - Satellites positions (in ECEF system: [X; Y; Z;] -
%                         one column per satellite) from Least Squares
%       transmitTime    - Transmitting time of channels to be processed 
%                         corresponding to measurement time 
%       localTime       - Local time(in GPST) at measurement time,
%                         uncorrected by Least Square's clock bias
%       settings        - Receiver settings 
%       satElev         - Satellites elevation angles (degrees) 
%                         from Least Squares
%       LS_clkBias      - Receiver clock bias estimated by Least Squares
%       satClkCorr      - Satellite clock bias

%
%   Outputs:
%       DPE_latitude    - Estimated receiver latitude from Corr-DPE
%       DPE_longitude   - Estimated receiver longitude from Corr-DPE
%       DPE_height      - Estimated receiver height from Corr-DPE
%       DPE_clkBias     - Estimated receiver clock bias from Corr-DPE

% -------------------------------------------------------------------------
% ------------------------ CorrDPE_ARS v1.0 --------------------------------
%   Copyright (C) 2025 Sergio Vicenzo
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
%
% Last Updated: 08 Dec 2025
% 


% === Initialize parameters ===============================================
dtr                     = pi/180;
trop                    = zeros(length(activeChnList),1);
chip_spacings           = [-(flip(settings.chipspacing_dpe_precalc:...
                            settings.chipspacing_dpe_precalc:settings.dllCorrelatorSpacing)),0,...
                            settings.chipspacing_dpe_precalc:...
                            settings.chipspacing_dpe_precalc:settings.dllCorrelatorSpacing];
precalc_correlations    = ...
                          zeros(length(chip_spacings)...
                          *length(activeChnList),3);

% === Make folder to store correlograms ===================================
if ~exist([settings.outfile_root '\Correlogram\'], 'dir') && ...
        settings.DPE_plotCorrelogram == 1
    mkdir([settings.outfile_root '\Correlogram\']);
end

if (settings.useTropCorr == 1)
    for i = 1:length(activeChnList)
     % === Estimate tropospheric error for each satellite =================
     trop(i) = tropo(sin(satElev(activeChnList(i)) * dtr), ...
         0.0, 1013.0, 293.0, 50.0, 0.0, 0.0, 0.0);
    end
end

% === Get multicorrelator values from tracking ============================
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

    % Generate indexes to store the multicorrelators in an array
    count2=3;
    count=1+(j-1)*length(chip_spacings);

 % === Pre-calculate the correlations =================================
 for closestIndex=nonCohInt
    % Loop repeats depending on the non-coherent integration time
    % chosen

    for ii=1:length(chip_spacings)
    
      % === Store the correlations and its corresponding pseudorange =====
      precalc_correlations(count,1) = ...
          trackResults(activeChnList(j)).PRN;
      precalc_correlations(count,2) = ...
          (localTime-transmitTime(j))*settings.c ...
          + (chip_spacings(ii)/settings.codeFreqBasis) * settings.c;
      precalc_correlations(count,count2) = ...
          trackResults(activeChnList(j)).MultiCorr_IQ(ii,closestIndex);

      count=count+1;
    end
      count=1+(j-1)*length(chip_spacings);
      count2=count2+1;
 end 
 % End for pre-calculating the corr values for a single satellite

end % End for pre-calculating the corr values for all satellites

% Non-coherent integration of each satellite's autocorrelation function
precalc_correlations(:,3) = sum(precalc_correlations(1:end,3:end),2);


%% === Initialize Accelerated Random Search (ARS) =========================
% The following code is adapted from the DPE simulator by Pau Closas and
% Tang Shuo in https://github.com/Shuo-Tang/direct_position_estimation
Niter = 2000;
navSolutions.Niter = Niter;
gamma_est = zeros(3,Niter+1);
contraction =2;          % contraction parameter
navSolutions.contraction = contraction;          % contraction parameter
dmax = 100;
navSolutions.dmax = dmax;
dmin = 1;
navSolutions.dmin = dmin;
dmax_clk=100;%1/settings.samplingFreq/10;
navSolutions.dmax_clk=dmax_clk;
dmin_clk=1;%1/settings.samplingFreq/100;
navSolutions.dmin_clk=dmin_clk;
d = dmax;
d_clk = dmax_clk;

amp_est = zeros(1,Niter+1);
EstRxClkBias =zeros(1,Niter+1);

% Corr-DPE estimates at each ARS iteration
% Position estimates
gamma_est(:,1) = [navSolutions.X(currMeasNr);...
  navSolutions.Y(currMeasNr);  navSolutions.Z(currMeasNr)]; 

% Clock bias estimates
EstRxClkBias(:,1)=LS_clkBias;

% Corr-DPE correlation values
corriter_prn_B_SAT_Sum=zeros(1,Niter+1);

% Initialize ARS iteration with LS PVT
for j=1:length(activeChnList)
  
    candia_xyz = gamma_est(1:3,1);
    
    candia_xyz_pseudoranges =  ...
        sqrt(sum((satPositions(:, j) - candia_xyz).^2)) ...
        + trop(j) - satClkCorr(j)* settings.c;
    clockk =  EstRxClkBias(:,1);
    
    candia_xyz_pseudoranges_clkBias...
        = candia_xyz_pseudoranges+clockk;

    prn_index = find(precalc_correlations(:,1)==...
            trackResults(activeChnList(j)).PRN);
               
     corriter_prn = interp1(precalc_correlations(prn_index,2),...
              precalc_correlations(prn_index,3),...
              candia_xyz_pseudoranges_clkBias,'linear');
     corriter_prn(isnan(corriter_prn)) = 0;
     corriter_prn_B_SAT_Sum(:,1) = ...
         corriter_prn+corriter_prn_B_SAT_Sum(:,1);

end

     J_ant = (corriter_prn_B_SAT_Sum(:,1));
     amp_est(:,1) = J_ant;
 
for it = 1:Niter        %%% ARS algorithm iterations
    
    % draw a random movement
    rand_point = gamma_est(:,it) + d*(2*rand(3,1)-1); % Random position
    rand_clk = EstRxClkBias(:,it)+d_clk*(2*rand-1); % Random clock bias

     for j=1:length(activeChnList)
      
      candia_xyz = rand_point;

        candia_xyz_pseudoranges =  ...
            sqrt(sum((satPositions(:, j) - candia_xyz).^2)) ...
            + trop(j) - satClkCorr(j)* settings.c;

        clockk =  rand_clk;
      
            candia_xyz_pseudoranges_clkBias...
                = candia_xyz_pseudoranges+clockk;
                   
        prn_index = find(precalc_correlations(:,1)==...
            trackResults(activeChnList(j)).PRN);
               
    
     corriter_prn=interp1(precalc_correlations(prn_index,2),...
              precalc_correlations(prn_index,3),...
              candia_xyz_pseudoranges_clkBias,'linear');
     corriter_prn(isnan(corriter_prn))=0;
     corriter_prn_B_SAT_Sum(:,it+1)=...
         corriter_prn+corriter_prn_B_SAT_Sum(:,it+1);
    

    end

     J = (corriter_prn_B_SAT_Sum(:,it+1));
    
    % select or discard point
    if J > J_ant 
        % If PVT bears higher correlation value than previous PVT guess, 
        % select this PVT as latest PVT estimate
        gamma_est(:,it+1) = rand_point;
        EstRxClkBias(:,it+1)=rand_clk;
        amp_est(:, it+1) =corriter_prn_B_SAT_Sum(:,it+1);
        J_ant = J;
        d = dmax;
        d_clk=dmax_clk;
    else
        % If PVT bears lower correlation value than previous PVT guess, 
        % contract ARS search space and use prev epoch's PVT estimate as
        % latest estimate
        gamma_est(:,it+1) = gamma_est(:,it);
        EstRxClkBias(:,it+1)=EstRxClkBias(:,it);
        amp_est(:,it+1) = amp_est(:,it);
        d = d/contraction;
        d_clk=d_clk/contraction;
    end
    if d < dmin
        d = dmax;
    end
    
    if d_clk < dmin_clk
        d_clk =dmax_clk;
    end
end

%% === Record DPE estimates ===============================================
[navSolutions.DPE_latitude(currMeasNr), ...
 navSolutions.DPE_longitude(currMeasNr), ...
 navSolutions.DPE_height(currMeasNr)] = cart2geo(...
                                    gamma_est(1,it+1), ...
                                    gamma_est(2,it+1), ...
                                    gamma_est(3,it+1), ...
                                    5);
navSolutions.DPE_estimate(currMeasNr, 1:5) =  ...
    [navSolutions.DPE_latitude(currMeasNr)...
    navSolutions.DPE_longitude(currMeasNr)...
    navSolutions.DPE_height(currMeasNr)...
    J_ant EstRxClkBias(:,it+1)] ;


navSolutions.DPE_clkBias(currMeasNr) = ...
    EstRxClkBias(:,it+1);

end


