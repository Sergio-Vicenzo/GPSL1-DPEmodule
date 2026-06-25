function [navSolutions]=DPE_module_PF7...
    (currMeasNr,navSolutions,activeChnList,trackResults,currMeasSample,...
    satPositions,~,localTime,settings,satElev,fid,...
    LS_clkBias,satClkCorr,satVelocity,codePhase,satClkCorrRat,~)

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
% [navSolutions]=DPE_module_PF7...
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
%       settings        - Receiver settings 
%       satElev         - Satellites elevation angles (degrees) 
%                         from Least Squares
%       fid             - Integer file identifier for raw IF data
%       LS_clkBias      - Receiver clock bias estimated by Least Squares
%       satClkCorr      - Satellite clock bias
%       satVelocity     - Satellites positions (in ECEF system: [X; Y; Z;] -
%                         one column per satellite) 
%       codePhase       - Estimated code phase of each signal from tracking
%       satClkCorrRat   - Satellite clock drift
%
%   Outputs:
%       DPE_estimate    - All of BDPE PVT estimates, including positions
%                           and velocity
%       DPE_latitude    - Estimated receiver latitude from DPE
%       DPE_longitude   - Estimated receiver longitude from DPE
%       DPE_height      - Estimated receiver height from DPE
%       DPE_clkBias     - Estimated receiver clock bias from DPE

% -------------------------------------------------------------------------
% ------------------------ BDPE_module v1.0 --------------------------------
%   Copyright (C) 2026 Sergio Vicenzo
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

% Last Updated: 11 June 2025

addpath include common

% === Initialize parameters ===============================================
dtr                     = pi/180;
m2lat                   = 1/110734;
m2lon                   = 1/103043;
GT_ECEF                 = llh2xyz(settings.gt_llh.*[pi/180 pi/180 1]);
trop                    = zeros(length(activeChnList),1);
iono                    = zeros(length(activeChnList),1);
caCode                  = zeros(length(activeChnList),settings.codeLength);
closestIndex            = zeros(length(activeChnList),1);
bit_idx                 = zeros(length(activeChnList),1);
navData_sign1           = zeros(length(activeChnList),1);
navData_sign2           = zeros(length(activeChnList),1);
numberOfparticles       = 10000;%20164;
record_correlogram_pos  = zeros(numberOfparticles,1);
record_correlogram_vel  = zeros(numberOfparticles,1);

% === Make folder to store correlograms ===================================
if ~exist([settings.outfile_root '\Correlogram\'], 'dir') && ...
        settings.DPE_plotCorrelogram == 1
    mkdir([settings.outfile_root '\Correlogram\']);
end

% === Generate Clock Bias search space ====================================
candia_xyz_LS = [navSolutions.X(currMeasNr),...
        navSolutions.Y(currMeasNr),navSolutions.Z(currMeasNr)];

% Generate C/A code for each satellite
for i = 1:length(activeChnList)
        caCode(i,:) = generateCAcode(trackResults(activeChnList(i)).PRN);

        [~,closestIndex(i)] = min(abs(currMeasSample-...
            trackResults(activeChnList(i)).absoluteSample));

        if currMeasSample < ...
              trackResults(activeChnList(i)).absoluteSample(closestIndex(i))
        
          closestIndex(i)=closestIndex(i)-1;
        end
end

if currMeasNr == 1

    % Least Square's clock bias estimates are found to be too inaccurate
    % for DPE. Thus, for the first epoch, DPE will first get a rought
    % estimate of the clock bias!

    clock_bias              = zeros(length(activeChnList),1);
    rough_codePhase         = zeros(length(activeChnList),1);
    carrFreq_real           = zeros(length(activeChnList),1);
    carrFreq                = zeros(length(activeChnList),1);
    clkDrift1               = zeros(length(activeChnList),1);
    clkDrift2               = zeros(length(activeChnList),1);

    for j = 1:length(activeChnList)

          if (settings.useTropCorr == 1)
            trop(j) = tropo(sin(satElev(activeChnList(j)) * dtr), ...
                 0.0, 1013.0, 293.0, 50.0, 0.0, 0.0, 0.0);
             navSolutions.trop(currMeasNr,j)=trop(j);
          end

          if (settings.useIonoCorr == 1)
              iono(j) = ionocorr(localTime,satPositions(:,j),...
                GT_ECEF,settings);
            navSolutions.iono(currMeasNr,j)=iono(j);
          end
       
        % Generate indexes to store the pre-calculated correlations in an
        % array
        rough_codePhase(j)=round(mod(codePhase(activeChnList(j))/1023,1)*1023,1);

        count3=1;

        clkBias_est=zeros(length(-100000:200000),2);
        for clkBias_space=-100000:200000
             candia_xyz_pseudorange =  ...
                    sqrt(sum(((satPositions(:,j)'...
                    - GT_ECEF).^2))) +...
                    trop(j) - satClkCorr(j)* settings.c + iono(j) + ...
                     LS_clkBias+clkBias_space;
        
               clkBias_est(count3,1)=...
                   settings.codeLength-mod(candia_xyz_pseudorange...
                        /(settings.c*1e-3),1)*settings.codeLength;
               clkBias_est(count3,2)=...
                   LS_clkBias+clkBias_space;
        
               count3=count3+1;
        end
  
        [~,anIndex]=min(abs(repmat(rough_codePhase(j),...
            length(-100000:200000),1)-clkBias_est(:,1)));
        
        clock_bias(j) = ...
            clkBias_est(anIndex,2);

        candia_vel= [0 0 0];

        carrFreq_real(j)= trackResults(activeChnList(j)).carrFreq(closestIndex(j));

         % Check for carrier inversion
         carrFreq(j) = (-1575.42e6/settings.c)*...
            ((satVelocity(:,j)-candia_vel')'*...
            ((satPositions(:,j)-GT_ECEF')/norm(satPositions(:,j)-GT_ECEF')) -...
            settings.c*satClkCorrRat(j) ) ;

         % With carrier inversion
         clkDrift1(j)=(carrFreq(j)+carrFreq_real(j))*settings.c/1575.42e6 ;

         % No carrier inversion
         clkDrift2(j)=(carrFreq(j)-carrFreq_real(j))*settings.c/1575.42e6 ;
              
    end

    rough_clkBias = mean(clock_bias(~isoutlier(clock_bias))); % in meters
    % std_clkBias = std(clock_bias);

    if var(clkDrift1) < var(clkDrift2) 
        % Carrier inversion
        navSolutions.carrInv = -1; %LABSAT
        rough_clkDrift = mean(clkDrift1); % in meters/second
    else 
        % No Carrier inversion
        navSolutions.carrInv = 1; %NSL STEREO
        rough_clkDrift = mean(clkDrift2); % in meters/second
    end

    navSolutions.rough_clkBias=rough_clkBias;
    navSolutions.rough_clkDrift=rough_clkDrift;

else % Get tropospheric and ionospheric corrections

    % DPE_pureClkBias=navSolutions.DPE_pureClkBias;

    if (settings.useTropCorr == 1)
        for i = 1:length(activeChnList)
         % === Estimate tropospheric error for each satellite =============
         trop(i) = tropo(sin(satElev(activeChnList(i)) * dtr), ...
             0.0, 1013.0, 293.0, 50.0, 0.0, 0.0, 0.0);
         navSolutions.trop(currMeasNr,i)=trop(i);
        end
    
    end

    if (settings.useIonoCorr == 1)
        for i = 1:length(activeChnList)
            %--- Calculate tropospheric correction --------------------
            iono(i) = ionocorr(localTime,satPositions(:,i),...
                GT_ECEF,settings);
            navSolutions.iono(currMeasNr,i)=iono(i);
        end
    
    end
    
end

if currMeasNr==1
    %% Initialize PVT particles and their corresponding weights
    posLat_particles    = navSolutions.latitude(currMeasNr)+...
                            m2lat.*normrnd(0,50,numberOfparticles,1);

    posLong_particles   = navSolutions.longitude(currMeasNr)+...
                            m2lon.*normrnd(0,50,numberOfparticles,1);

    posHeight_particles = navSolutions.height(currMeasNr)+...
                                   normrnd(0,50,numberOfparticles,1);
    % 
    % posLat_particles    = settings.gt_llh(1)+...
    %                         m2lat.*normrnd(0,20,numberOfparticles,1);
    % 
    % posLong_particles   = settings.gt_llh(2)+...
    %                         m2lon.*normrnd(0,20,numberOfparticles,1);
    % 
    % posHeight_particles = settings.gt_llh(3)+...
    %                                normrnd(0,20,numberOfparticles,1);

    x = posLat_particles./180 .* pi;
    y = posLong_particles./180 .* pi;

    candia_xyz  = llh2xyz([x,y,posHeight_particles]); 

    posX_particles = candia_xyz(:,1);
    posY_particles = candia_xyz(:,2);
    posZ_particles = candia_xyz(:,3);

    % clkBias_particles   = rough_clkBias+...
    %                         normrnd(0,50,numberOfparticles,1); % in m
    clkBias_particles   = normrnd(0,50,numberOfparticles,1); % in m
    Xvelocity_particles = normrnd(0,10,numberOfparticles,1);
    Yvelocity_particles = normrnd(0,10,numberOfparticles,1);
    Zvelocity_particles = normrnd(0,10,numberOfparticles,1);
    clkDrift_particles  = rough_clkDrift+...
                            normrnd(0,10,numberOfparticles,1); % in m/s
    
    % Save the particles for next positioning epoch
    navSolutions.posX_particles         = posX_particles;
    navSolutions.posY_particles         = posY_particles;
    navSolutions.posZ_particles         = posZ_particles;
    navSolutions.clkBias_particles      = clkBias_particles;
    navSolutions.Xvelocity_particles    = Xvelocity_particles;
    navSolutions.Yvelocity_particles    = Yvelocity_particles;
    navSolutions.Zvelocity_particles    = Zvelocity_particles;
    navSolutions.clkDrift_particles     = clkDrift_particles;

    weight_pos  = ones(1,numberOfparticles)/numberOfparticles;
    weight_vel  = ones(1,numberOfparticles)/numberOfparticles;

else
    %% Propagate the particles from the previous positioning epoch

    % % Establish process noise, Q
    % process_noise(1:3)  = ones(3,1)*5;%((settings.navSolPeriod/1000)^3./3).*(10^-2);
    % process_noise(4)    = 5;%1e-1;
    % process_noise(5:7)  = ones(3,1)*1;%((settings.navSolPeriod/1000)).*(10^-2);
    % process_noise(8)    = 1;%1e-2;

    % Establish process noise, Q
    process_noise(1:3)  = ones(3,1)*2;%((settings.navSolPeriod/1000)^3./3).*(10^-2);
    process_noise(4)    = 2;%1e-1;
    process_noise(5:7)  = ones(3,1)*0.5;%((settings.navSolPeriod/1000)).*(10^-2);
    process_noise(8)    = 0.5;%1e-2;

    

    for i = 1:numberOfparticles
    
        posX_particles(i,1)      = normrnd(navSolutions.posX_particles(i),...
                                    process_noise(1),1,1) + ...
                                    (settings.navSolPeriod/1000)*...
                                    navSolutions.DPE_estimate(currMeasNr-1,5);

        posY_particles(i,1)      = normrnd(navSolutions.posY_particles(i),...
                                    process_noise(2),1,1) + ...
                                    (settings.navSolPeriod/1000)*...
                                    navSolutions.DPE_estimate(currMeasNr-1,6);

        posZ_particles(i,1)      = normrnd(navSolutions.posZ_particles(i),...
                                    process_noise(3),1,1) + ...
                                    (settings.navSolPeriod/1000)*...
                                    navSolutions.DPE_estimate(currMeasNr-1,7);

        % clkBias_particles(i)   = normrnd(navSolutions.clkBias_particles(i),...
        %                             process_noise(4),1,1) + ...
        %                             (settings.navSolPeriod/1000)*...
        %                             navSolutions.DPE_estimate(currMeasNr-1,8);

        clkBias_particles(i,1)   = normrnd(navSolutions.clkBias_particles(i),...
                                    process_noise(4),1,1) ;

        Xvelocity_particles(i,1) = normrnd(navSolutions.Xvelocity_particles(i),...
                                    process_noise(5),1,1);

        Yvelocity_particles(i,1) = normrnd(navSolutions.Yvelocity_particles(i),...
                                    process_noise(6),1,1);
        
        Zvelocity_particles(i,1) = normrnd(navSolutions.Zvelocity_particles(i),...
                                    process_noise(7),1,1);

        clkDrift_particles(i,1)  = normrnd(navSolutions.clkDrift_particles(i),...
                                    process_noise(8),1,1);

    end

    % Save the particles for next positioning epoch
    navSolutions.posX_particles         = posX_particles;
    navSolutions.posY_particles         = posY_particles;
    navSolutions.posZ_particles         = posZ_particles;
    navSolutions.clkBias_particles      = clkBias_particles;
    navSolutions.Xvelocity_particles    = Xvelocity_particles;
    navSolutions.Yvelocity_particles    = Yvelocity_particles;
    navSolutions.Zvelocity_particles    = Zvelocity_particles;
    navSolutions.clkDrift_particles     = clkDrift_particles;

    % Use the weights from previous positioning epoch
    weight_pos = navSolutions.weight_pos(currMeasNr-1,:);
    weight_vel = navSolutions.weight_vel(currMeasNr-1,:);

end

% Precalc correlations across code delay and Doppler freq
code_spacings                   = zeros(length( [-300:1:300]));
Dopp_spacings                   = zeros(length([-50:1:50]));
precalc_correlations_Doppler    = zeros(length(Dopp_spacings)...
                                    *length(activeChnList),3);
precalc_correlations_codeDelay  = zeros(length(code_spacings)...
                                    *length(activeChnList),3);


if currMeasNr == 1
    clk_Bias    = rough_clkBias;
    clk_Drift   = rough_clkDrift;
    prev_PosEst = GT_ECEF;
    prev_VelEst = [0 0 0];
else
    clk_Bias    = navSolutions.DPE_estimate(currMeasNr-1,4)...
                    +(settings.navSolPeriod/1000)*...
                    navSolutions.DPE_estimate(currMeasNr-1,8);
    clk_Drift   = navSolutions.DPE_estimate(currMeasNr-1,8);
    prev_PosEst = navSolutions.DPE_estimate(currMeasNr-1,1:3)+ ...
                            (settings.navSolPeriod/1000)*...
                            navSolutions.DPE_estimate(currMeasNr-1,5:7);
    prev_VelEst = navSolutions.DPE_estimate(currMeasNr-1,5:7);
end

count5=1;
count6=1;

for i = 1:length(activeChnList)

    carrFreq = navSolutions.carrInv*(-1575.42e6/settings.c)*...
            ((satVelocity(:,i)-prev_VelEst')'*...
            ((satPositions(:,i)-prev_PosEst')/norm(satPositions(:,i)-prev_PosEst'))+...
            clk_Drift -...
            settings.c*satClkCorrRat(i) ) ;

    codeFreq = settings.codeFreqBasis ...
            + navSolutions.carrInv*(carrFreq*settings.codeFreqBasis / 1575.42e6);
    
    codePhaseStep = ...
            codeFreq...
            / settings.samplingFreq;

    blksize = round((settings.codeLength*settings.DPE_cohInt) ...
             / codePhaseStep);

    bit_idx(i) = ...
            trackResults(activeChnList(i)).absoluteSample(closestIndex(i))...
            + blksize - currMeasSample;      % 391238 samples (15.05 ms)

    if bit_idx(i)<=1
         bit_idx(i) = ...
            trackResults(activeChnList(i)).absoluteSample(closestIndex(i)+1)...
            - currMeasSample;      % 391238 samples (15.05 ms)
    end

    if bit_idx(i)==0
        bit_idx(i) = 2;
        % navData_sign2(j)=navData_sign1(j);
    end

    if settings.DPE_cohInt == 1
        bit_idx(i)=2;
        navData_sign1(i)=1;
        navData_sign2(i)=navData_sign1(i);
    end

    navData_sign1(i) = sign(trackResults(activeChnList(i)).I_P(closestIndex(i)));
    navData_sign2(i) = sign(trackResults(activeChnList(i)).I_P(closestIndex(i)+1));


      % === Record sample number (based on 8-bit and 16-bit samples) ====
    if strcmp(settings.dataType,'int16') 
        fseek(fid, ...
            (settings.fileType*...
            (currMeasSample)*2),'bof');
    else
        fseek(fid, ...
            settings.fileType*...
            (currMeasSample),'bof');
    end

     % === Read in the appropriate number of samples to process =======
    [rawSignal, ~] = fread(fid, ...
                   settings.fileType*blksize, settings.dataType);
        
    rawSignal = rawSignal';
    
    if (settings.fileType==2)
       rawSignal1=rawSignal(1:2:end);
       rawSignal2=rawSignal(2:2:end);
       rawSignal = rawSignal1 + 1i .* rawSignal2;%transpose vector
    end

    code_spacings(:,i)= (sqrt(sum(((satPositions(:,i)'...
                        - prev_PosEst).^2))) +...
                        trop(i) - satClkCorr(i)* settings.c + iono(i) +...
                        clk_Bias) +  [-300:1:300];

    Dopp_spacings(:,i)=navSolutions.carrInv*(-1575.42e6/settings.c)*...
                        ((satVelocity(:,i)-prev_VelEst')'*...
                        ((satPositions(:,i)-prev_PosEst')/...
                        norm(satPositions(:,i)-prev_PosEst'))+...
                        clk_Drift-...
                        settings.c*satClkCorrRat(i)) +...
                        [-50:1:50];

    candia_xyz_pseudorange =  ...
                        sqrt(sum(((satPositions(:,i)'...
                        - prev_PosEst).^2))) +...
                        trop(i) - satClkCorr(i)* settings.c + iono(i) + ...
                        clk_Bias;

   Frac_codeDelay =...
               settings.codeLength...
                    -mod(candia_xyz_pseudorange/(settings.c*1e-3),1)*1023;

   delay_index = ...
                Frac_codeDelay: ...
                codePhaseStep : ...
                ((blksize-1)*codePhaseStep + ...
                Frac_codeDelay);       

    caCode1 = [ caCode(i,end) repmat(caCode(i,:),1,1)...
                repmat(caCode(i,:),1,settings.DPE_cohInt) ...
                repmat(caCode(i,:),1,1) caCode(i,1)];
    tcodee = ceil(delay_index)+1;

    s = caCode1(tcodee);
    
    % Precalc Doppler freq first
    for j = 1:length(Dopp_spacings(:,i))

        carrFreq = Dopp_spacings(j,i);

        % === Get the argument to sin/cos functions =======================
        time    = (0:blksize) ./ settings.samplingFreq;

        trigarg = ...
            ((carrFreq...
            * 2.0 * pi) .* time);

        % === Compute the carrier signal ==================================
        carrsig = exp(1i .* trigarg(1:blksize));

        % === Mix the carrier replica to baseband =========================
        qBasebandSignal1 = real(carrsig(1:bit_idx(i)-1) .* rawSignal(1:bit_idx(i)-1));
        iBasebandSignal1 = imag(carrsig(1:bit_idx(i)-1) .* rawSignal(1:bit_idx(i)-1));

        qBasebandSignal2 = real(carrsig(bit_idx(i):end) .* rawSignal(bit_idx(i):end));
        iBasebandSignal2 = imag(carrsig(bit_idx(i):end) .* rawSignal(bit_idx(i):end));

    
        I1 = sum(navData_sign1(i).*s(1:bit_idx(i)-1) .* iBasebandSignal1);
        Q1 = sum(navData_sign1(i).*s(1:bit_idx(i)-1) .* qBasebandSignal1);

        I2 = sum(navData_sign2(i).*s(bit_idx(i):end) .* iBasebandSignal2);
        Q2 = sum(navData_sign2(i).*s(bit_idx(i):end) .* qBasebandSignal2);

        % correlogram_single=(I1 + I2)^2 + (Q1 + Q2)^2;
         precalc_correlations_Doppler(count5,1) = ...
              trackResults(activeChnList(i)).PRN;
         precalc_correlations_Doppler(count5,2) = ...
              Dopp_spacings(j,i);
         precalc_correlations_Doppler(count5,3)=(I1 + I2)^2 + (Q1 + Q2)^2;

         count5=count5+1;

    end

    carrFreq = navSolutions.carrInv*(-1575.42e6/settings.c)*...
            ((satVelocity(:,i)-prev_VelEst')'*...
            ((satPositions(:,i)-prev_PosEst')/norm(satPositions(:,i)-prev_PosEst'))+...
            clk_Drift -...
            settings.c*satClkCorrRat(i) ) ;

    % === Get the argument to sin/cos functions =======================
    time    = (0:blksize) ./ settings.samplingFreq;

    trigarg = ...
        ((carrFreq...
        * 2.0 * pi) .* time);

    % === Compute the carrier signal ==================================
    carrsig = exp(1i .* trigarg(1:blksize));

    % === Mix the carrier replica to baseband =========================
    qBasebandSignal1 = real(carrsig(1:bit_idx(i)-1) .* rawSignal(1:bit_idx(i)-1));
    iBasebandSignal1 = imag(carrsig(1:bit_idx(i)-1) .* rawSignal(1:bit_idx(i)-1));

    qBasebandSignal2 = real(carrsig(bit_idx(i):end) .* rawSignal(bit_idx(i):end));
    iBasebandSignal2 = imag(carrsig(bit_idx(i):end) .* rawSignal(bit_idx(i):end));

    for j = 1:length(code_spacings(:,i))

        Frac_codeDelay=...
           settings.codeLength...
                -mod(code_spacings(j,i)/(settings.c*1e-3),1)*1023;

         delay_index = ...
                Frac_codeDelay: ...
                codePhaseStep : ...
                ((blksize-1)*codePhaseStep + ...
                Frac_codeDelay);       
        
        tcodee = ceil(delay_index)+1;

        s = caCode1(tcodee);

        I1 = sum(navData_sign1(i).*s(1:bit_idx(i)-1) .* iBasebandSignal1);
        Q1 = sum(navData_sign1(i).*s(1:bit_idx(i)-1) .* qBasebandSignal1);

        I2 = sum(navData_sign2(i).*s(bit_idx(i):end) .* iBasebandSignal2);
        Q2 = sum(navData_sign2(i).*s(bit_idx(i):end) .* qBasebandSignal2);

         precalc_correlations_codeDelay(count6,1) = ...
              trackResults(activeChnList(i)).PRN;
         precalc_correlations_codeDelay(count6,2) = ...
              code_spacings(j,i);
         precalc_correlations_codeDelay(count6,3)=(I1 + I2)^2 + (Q1 + Q2)^2;

         count6=count6+1;

    end
    prn_index = find(precalc_correlations_codeDelay(:,1)==...
            trackResults(activeChnList(i)).PRN);

    [~,pseudo_index]=max(precalc_correlations_codeDelay(prn_index,3));

    temp_precalc_correlations=precalc_correlations_codeDelay(prn_index,:);

    navSolutions.estCodeDelay(currMeasNr,activeChnList(i))=...
        temp_precalc_correlations(pseudo_index,2)-clk_Bias...
            +satClkCorr(i)* settings.c ;

end

for i = 1:numberOfparticles

    % === Pre-allocate space for correlogram ==============================
    correlogram=0;

    % Get particle PVTs
    candia_xyz = [posX_particles(i) posY_particles(i) posZ_particles(i)];


    for j=1:length(activeChnList)

        candia_xyz_pseudorange =  ...
            sqrt(sum(((satPositions(:,j)'...
            - candia_xyz).^2))) +...
            trop(j) - satClkCorr(j)* settings.c + iono(j) ...
            + clk_Bias ...
            + clkBias_particles(i);

        prn_index = find(precalc_correlations_codeDelay(:,1)==...
            trackResults(activeChnList(j)).PRN);
    
            % === Obtain correlogram from single satellite ====================
        correlogram_single = ...
        interp1(precalc_correlations_codeDelay(prn_index,2),...
          precalc_correlations_codeDelay(prn_index,3),candia_xyz_pseudorange,'linear');

        if isnan(correlogram_single)
            correlogram_single=0;
        end
         
        correlogram=correlogram+correlogram_single;

    end

    record_correlogram_pos(i)=correlogram;

    weight_exp=0;

    for exponentialSum_idx=0:20
        weight_exp=weight_exp+...
            (correlogram^exponentialSum_idx)/factorial(exponentialSum_idx);
    end

    weight_pos(i)=(weight_pos(i))*weight_exp;
    
end

for i = 1:numberOfparticles

    % === Pre-allocate space for correlogram ==============================
    correlogram=0;

    % Get particle PVTs
    % candia_xyz = [posX_particles(i) posY_particles(i) posZ_particles(i)];
    candia_xyz = prev_PosEst;
    candia_vel = [Xvelocity_particles(i) ...
                        Yvelocity_particles(i) Zvelocity_particles(i)];

    for j=1:length(activeChnList)

        carrFreq = navSolutions.carrInv*(-1575.42e6/settings.c)*...
            ((satVelocity(:,j)-candia_vel')'*...
            ((satPositions(:,j)-candia_xyz')/norm(satPositions(:,j)-candia_xyz'))+...
            clkDrift_particles(i) -...
            settings.c*satClkCorrRat(j) ) ;

         prn_index = find(precalc_correlations_Doppler(:,1)==...
            trackResults(activeChnList(j)).PRN);
    
        % === Obtain correlogram from single satellite ====================
        correlogram_single = ...
            interp1(precalc_correlations_Doppler(prn_index,2),...
              precalc_correlations_Doppler(prn_index,3),carrFreq,'linear');

        if isnan(correlogram_single)
            correlogram_single=0;
        end
         
        correlogram=correlogram+correlogram_single;

    end

    record_correlogram_vel(i)=correlogram;

    % Approximate exponential weights
    weight_exp=0;

    for exponentialSum_idx=0:20
        weight_exp=weight_exp+...
            (correlogram^exponentialSum_idx)/factorial(exponentialSum_idx);
    end

    weight_vel(i)=(weight_vel(i))*weight_exp;
    
end


% Normalise weights
weight_pos=weight_pos./sum(weight_pos);
weight_vel=weight_vel./sum(weight_vel);

% Get PF-DPE estimate
DPE_X       = sum(weight_pos'.*posX_particles);
DPE_Y       = sum(weight_pos'.*posY_particles);
DPE_Z       = sum(weight_pos'.*posZ_particles);
DPE_clkBias = sum(weight_pos'.*clkBias_particles) + clk_Bias;

DPE_Xvel    = sum(weight_vel'.*Xvelocity_particles);
DPE_Yvel    = sum(weight_vel'.*Yvelocity_particles);
DPE_Zvel    = sum(weight_vel'.*Zvelocity_particles);
DPE_clkDrift= sum(weight_vel'.*clkDrift_particles);


% Calculate effective no. of particles
N_eff=1/sum(weight_pos.^2);

% Resample with multinominal sampling
if N_eff < (2/3)*numberOfparticles

    Q   = [0 cumsum(weight_pos)];           % cumulative sum
    u   = rand(1,numberOfparticles);               % N uniform(0,1]
    idx = zeros(1,numberOfparticles);
    for n = 1:numberOfparticles
        idx(n) = find(Q >= u(n), 1) - 1;  % linear search
    end
    
    % Save the particles for next positioning epoch
    navSolutions.posX_particles         = posX_particles(idx);
    navSolutions.posY_particles         = posY_particles(idx);
    navSolutions.posZ_particles         = posZ_particles(idx);
    navSolutions.clkBias_particles      = clkBias_particles(idx);
    
    weight_pos=ones(1,numberOfparticles)/numberOfparticles;
end

navSolutions.weight_pos(currMeasNr,:)=weight_pos;

% Calculate effective no. of particles
N_eff=1/sum(weight_vel.^2);

% Resample with multinominal sampling
if N_eff < (2/3)*numberOfparticles

    Q   = [0 cumsum(weight_vel)];           % cumulative sum
    u   = rand(1,numberOfparticles);               % N uniform(0,1]
    idx = zeros(1,numberOfparticles);
    for n = 1:numberOfparticles
        idx(n) = find(Q >= u(n), 1) - 1;  % linear search
    end
    
    % Save the particles for next positioning epoch
    navSolutions.Xvelocity_particles    = Xvelocity_particles(idx);
    navSolutions.Yvelocity_particles    = Yvelocity_particles(idx);
    navSolutions.Zvelocity_particles    = Zvelocity_particles(idx);
    navSolutions.clkDrift_particles     = clkDrift_particles(idx);
    
    weight_vel=ones(1,numberOfparticles)/numberOfparticles;
end

navSolutions.weight_vel(currMeasNr,:)=weight_vel;


% === Record DPE estimates ================================================
navSolutions.DPE_estimate(currMeasNr, 1:8) =  ...
    [DPE_X DPE_Y DPE_Z DPE_clkBias DPE_Xvel DPE_Yvel DPE_Zvel DPE_clkDrift];

DPE_latlongheight=xyz2llh([DPE_X DPE_Y DPE_Z]);

DPE_latlongheight(1:2)=DPE_latlongheight(1:2).*180./pi;

navSolutions.DPE_latitude(currMeasNr) =  ...
    DPE_latlongheight(1);

navSolutions.DPE_longitude(currMeasNr) =  ...
    DPE_latlongheight(2);

navSolutions.DPE_height(currMeasNr) =  ...
    DPE_latlongheight(3);

navSolutions.DPE_clkBias(currMeasNr) = ...
    DPE_clkBias;


end


