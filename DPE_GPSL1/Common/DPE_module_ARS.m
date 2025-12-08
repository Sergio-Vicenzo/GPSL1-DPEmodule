function [navSolutions]=DPE_module_ARS...
    (currMeasNr,navSolutions,activeChnList,trackResults,currMeasSample,...
    satPositions,~,localTime,settings,satElev,fid,...
    LS_clkBias,satClkCorr,satVelocity,codePhase,satClkCorrRat)

% DPE_module_ARS is a Direct Position Estimation (DPE) plug-in module that can 
% be integrated into existing two-step positioning MATLAB SDRs.
%
%
% [navSolutions]=DPE_module_ARS...
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

%
%   Outputs:
%       DPE_latitude    - Estimated receiver latitude from DPE
%       DPE_longitude   - Estimated receiver longitude from DPE
%       DPE_height      - Estimated receiver height from DPE
%       DPE_clkBias     - Estimated receiver clock bias from DPE

% -------------------------------------------------------------------------
% ------------------------ DPE_module_ARS v1.0 --------------------------------
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

% Least Square's clock bias estimates are found to be too inaccurate
% for DPE. Thus, for the first epoch, DPE will first get a rought
% estimate of the clock bias!

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

      if (settings.useTropCorr == 1)
          iono(j) = ionocorr(localTime,satPositions(:,j),...
            candia_xyz_LS,settings);
        navSolutions.iono(currMeasNr,j)=iono(j);
      end
   

    % Generate indexes to store the pre-calculated correlations in an
    % array
    rough_codePhase(j)=mod(codePhase(activeChnList(j))/1023,1)*1023;

    count3=1;

    clkBias_est=zeros(length(-100000:100000),2);
    for clkBias_space=-100000:100000
         candia_xyz_pseudorange =  ...
                sqrt(sum(((satPositions(:,j)'...
                - candia_xyz_LS).^2))) +...
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
        length(-100000:100000),1)-clkBias_est(:,1)));
    
    clock_bias(j) = ...
        clkBias_est(anIndex,2);

    candia_vel= [0 0 0];

    carrFreq_real(j)= trackResults(activeChnList(j)).carrFreq(closestIndex(j));


     % Check for carrier inversion
     carrFreq(j) = (-1575.42e6/settings.c)*...
        ((satVelocity(:,j)-candia_vel')'*...
        ((satPositions(:,j)-candia_xyz_LS')/norm(satPositions(:,j)-candia_xyz_LS')) -...
        settings.c*satClkCorrRat(j) ) ;

     % With carrier inversion
     clkDrift1(j)=(carrFreq(j)+carrFreq_real(j))*settings.c/1575.42e6 ;

     % No carrier inversion
     clkDrift2(j)=(carrFreq(j)-carrFreq_real(j))*settings.c/1575.42e6 ;
 
end

rough_clkBias = mean(clock_bias); % in meters

if var(clkDrift1) < var(clkDrift2) 
    % Carrier inversion
    navSolutions.carrInv = -1; %LABSAT
    rough_clkDrift = mean(clkDrift1); % in meters/second
else 
    % No Carrier inversion
    navSolutions.carrInv = 1; %NSL STEREO
    rough_clkDrift = mean(clkDrift2); % in meters/second
end    

navSolutions.rough_clkDrift(currMeasNr)= rough_clkDrift;
navSolutions.rough_clkBias(currMeasNr) = rough_clkBias;

Niter=2000;
navSolutions.Niter=Niter;
gamma_est=zeros(3,Niter+1);
contraction =2;          % contraction parameter
navSolutions.contraction =contraction;          % contraction parameter
dmax = 1000;
navSolutions.dmax = dmax;
dmin = 1;
navSolutions.dmin = dmin;
dmax_clk=1000;%std_clkBias*3;%1/settings.samplingFreq/10;
navSolutions.dmax_clk=dmax_clk;
dmin_clk=1;%1/settings.samplingFreq/100;
navSolutions.dmin_clk=dmin_clk;
dmax_vel = 1;
navSolutions.dmax_vel = dmax_vel;
dmin_vel = 0.001;
navSolutions.dmin_vel = dmin_vel;
dmax_clkDrift = 1;
navSolutions.dmax_clkDrift = dmax_clkDrift;
dmin_clkDrift = 0.001;
navSolutions.dmin_clkDrift = dmin_clkDrift;
d = dmax;
d_clk=dmax_clk;
d_vel = dmax_vel;
d_clkDrift=dmax_clkDrift;

amp_est = zeros(1,Niter+1);
EstRxClkBias =zeros(1,Niter+1);

gamma_est(:,1) = [navSolutions.X(currMeasNr);...
  navSolutions.Y(currMeasNr);  navSolutions.Z(currMeasNr)]; 
gamma_vel(:,1) = [0 0 0];
EstRxClkBias(:,1)=rough_clkBias;
ExtRxClkDrift(:,1)=rough_clkDrift;
corriter_prn_B_SAT_Sum=zeros(1,Niter+1);


candia_xyz = gamma_est(1:3,1)';
candia_vel = gamma_vel(1:3,1)';

for j=1:length(activeChnList)


        carrFreq = navSolutions.carrInv*(-1575.42e6/settings.c)*...
            ((satVelocity(:,j)-candia_vel')'*...
            ((satPositions(:,j)-candia_xyz')/norm(satPositions(:,j)-candia_xyz'))+...
            ExtRxClkDrift(:,1) -...
            settings.c*satClkCorrRat(j) ) ;

        codeFreq = settings.codeFreqBasis ...
            + navSolutions.carrInv*(carrFreq*settings.codeFreqBasis / 1575.42e6);
    
        codePhaseStep = ...
                codeFreq...
                / settings.samplingFreq;
    
        blksize = round((settings.codeLength*settings.DPE_cohInt) ...
                 / codePhaseStep);
    
        bit_idx(j) = ...
                trackResults(activeChnList(j)).absoluteSample(closestIndex(j))...
                + blksize - currMeasSample;      % 391238 samples (15.05 ms)
    
        if bit_idx(j)<=1
             bit_idx(j) = ...
                trackResults(activeChnList(j)).absoluteSample(closestIndex(j)+1)...
                - currMeasSample;      % 391238 samples (15.05 ms)
        end
    
        if bit_idx(j)==0
            bit_idx(j) = 2;
            % navData_sign2(j)=navData_sign1(j);
        end
    
        if settings.DPE_cohInt == 1
            bit_idx(j)=2;
            navData_sign1(j)=1;
            navData_sign2(j)=navData_sign1(j);
        end
    
        navData_sign1(j) = sign(trackResults(activeChnList(j)).I_P(closestIndex(j)));
        navData_sign2(j) = sign(trackResults(activeChnList(j)).I_P(closestIndex(j)+1));

        
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
        

        % === Get the argument to sin/cos functions =======================
        time    = (0:blksize) ./ settings.samplingFreq;

        trigarg = ...
            ((carrFreq...
            * 2.0 * pi) .* time);

        % === Compute the carrier signal ==================================
        carrsig = exp(1i .* trigarg(1:blksize));

        % === Mix the carrier replica to baseband =========================
        qBasebandSignal1 = real(carrsig(1:bit_idx(j)-1) .* rawSignal(1:bit_idx(j)-1));
        iBasebandSignal1 = imag(carrsig(1:bit_idx(j)-1) .* rawSignal(1:bit_idx(j)-1));

        qBasebandSignal2 = real(carrsig(bit_idx(j):end) .* rawSignal(bit_idx(j):end));
        iBasebandSignal2 = imag(carrsig(bit_idx(j):end) .* rawSignal(bit_idx(j):end));

        candia_xyz_pseudorange =  ...
            sqrt(sum(((satPositions(:,j)'...
            - candia_xyz).^2))) +...
            trop(j) - satClkCorr(j)* settings.c + iono(j) + ...
            EstRxClkBias(:,1);

         Frac_codeDelay=...
           settings.codeLength...
                -mod(candia_xyz_pseudorange/(settings.c*1e-3),1)*1023;
         % Code period of L1 C/A is 1e-3 seconds

         delay_index = ...
                Frac_codeDelay: ...
                codePhaseStep : ...
                ((blksize-1)*codePhaseStep + ...
                Frac_codeDelay);       
        
        caCode1 = [ caCode(j,end) repmat(caCode(j,:),1,1)...
            repmat(caCode(j,:),1,settings.DPE_cohInt) ...
            repmat(caCode(j,:),1,1) caCode(j,1)];
        tcodee = ceil(delay_index)+1024;

        s = caCode1(tcodee);

        I1 = sum(navData_sign1(j).*s(1:bit_idx(j)-1) .* iBasebandSignal1);
        Q1 = sum(navData_sign1(j).*s(1:bit_idx(j)-1) .* qBasebandSignal1);

        I2 = sum(navData_sign2(j).*s(bit_idx(j):end) .* iBasebandSignal2);
        Q2 = sum(navData_sign2(j).*s(bit_idx(j):end) .* qBasebandSignal2);
           
               
        correlogram_single=(I1 + I2)^2 + (Q1 + Q2)^2;
         
        corriter_prn_B_SAT_Sum(:,1)=corriter_prn_B_SAT_Sum(:,1)+correlogram_single;

end

     J_ant = (corriter_prn_B_SAT_Sum(:,1));
     amp_est(:,1) = J_ant;

 
for it = 1:Niter        %%% ARS algorithm iterations
    
    % draw a random movement
    rand_point = gamma_est(:,it) + d*(2*rand(3,1)-1);
    rand_vel = [0;0;0];
    rand_clk = EstRxClkBias(:,it)+d_clk*(2*rand-1);
    rand_clkDrift = rough_clkDrift;%ExtRxClkDrift(:,it)+d_clkDrift*(2*rand-1);

     for j=1:length(activeChnList)

         carrFreq = navSolutions.carrInv*(-1575.42e6/settings.c)*...
                    ((satVelocity(:,j)-candia_vel')'*...
                    ((satPositions(:,j)-candia_xyz')/norm(satPositions(:,j)-candia_xyz'))+...
                    ExtRxClkDrift(:,1) -...
                    settings.c*satClkCorrRat(j) ) ;

        codeFreq = settings.codeFreqBasis ...
            + navSolutions.carrInv*(carrFreq*settings.codeFreqBasis / 1575.42e6);
    
        codePhaseStep = ...
                codeFreq...
                / settings.samplingFreq;
    
        blksize = round((settings.codeLength*settings.DPE_cohInt) ...
                 / codePhaseStep);
    
        bit_idx(j) = ...
                trackResults(activeChnList(j)).absoluteSample(closestIndex(j))...
                + blksize - currMeasSample;      % 391238 samples (15.05 ms)
    
        if bit_idx(j)<=1
             bit_idx(j) = ...
                trackResults(activeChnList(j)).absoluteSample(closestIndex(j)+1)...
                - currMeasSample;      % 391238 samples (15.05 ms)
        end
    
        if bit_idx(j)==0
            bit_idx(j) = 2;
            % navData_sign2(j)=navData_sign1(j);
        end
    
        if settings.DPE_cohInt == 1
            bit_idx(j)=2;
            navData_sign1(j)=1;
            navData_sign2(j)=navData_sign1(j);
        end
    
        navData_sign1(j) = sign(trackResults(activeChnList(j)).I_P(closestIndex(j)));
        navData_sign2(j) = sign(trackResults(activeChnList(j)).I_P(closestIndex(j)+1));

        
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
        

        % === Get the argument to sin/cos functions =======================
        time    = (0:blksize) ./ settings.samplingFreq;

        trigarg = ...
            ((carrFreq...
            * 2.0 * pi) .* time);

        % === Compute the carrier signal ==================================
        carrsig = exp(1i .* trigarg(1:blksize));

        % === Mix the carrier replica to baseband =========================
        qBasebandSignal1 = real(carrsig(1:bit_idx(j)-1) .* rawSignal(1:bit_idx(j)-1));
        iBasebandSignal1 = imag(carrsig(1:bit_idx(j)-1) .* rawSignal(1:bit_idx(j)-1));

        qBasebandSignal2 = real(carrsig(bit_idx(j):end) .* rawSignal(bit_idx(j):end));
        iBasebandSignal2 = imag(carrsig(bit_idx(j):end) .* rawSignal(bit_idx(j):end));

        candia_xyz_pseudorange =  ...
            sqrt(sum(((satPositions(:,j)...
            - rand_point).^2))) +...
            trop(j) - satClkCorr(j)* settings.c + iono(j) + ...
            rand_clk;

         Frac_codeDelay=...
           settings.codeLength...
                -mod(candia_xyz_pseudorange/(settings.c*1e-3),1)*1023;
         % Code period of L1 C/A is 1e-3 seconds

         delay_index = ...
                Frac_codeDelay: ...
                codePhaseStep : ...
                ((blksize-1)*codePhaseStep + ...
                Frac_codeDelay);       
        
        caCode1 = [ caCode(j,end) repmat(caCode(j,:),1,1)...
            repmat(caCode(j,:),1,settings.DPE_cohInt) ...
            repmat(caCode(j,:),1,1) caCode(j,1)];
        tcodee = ceil(delay_index)+1024;

        s = caCode1(tcodee);

        I1 = sum(navData_sign1(j).*s(1:bit_idx(j)-1) .* iBasebandSignal1);
        Q1 = sum(navData_sign1(j).*s(1:bit_idx(j)-1) .* qBasebandSignal1);

        I2 = sum(navData_sign2(j).*s(bit_idx(j):end) .* iBasebandSignal2);
        Q2 = sum(navData_sign2(j).*s(bit_idx(j):end) .* qBasebandSignal2);
           
               
        correlogram_single=(I1 + I2)^2 + (Q1 + Q2)^2;
         
        corriter_prn_B_SAT_Sum(:,it+1)=corriter_prn_B_SAT_Sum(:,it+1)+correlogram_single;

    end

     J = (corriter_prn_B_SAT_Sum(:,it+1));
    
    % select or discard point
    if J > J_ant
        gamma_est(:,it+1) = rand_point;
        gamma_vel(:,it+1) = rand_vel;

        ExtRxClkDrift(:,it+1)=rand_clkDrift;
        EstRxClkBias(:,it+1)=rand_clk;
        amp_est(:, it+1) =corriter_prn_B_SAT_Sum(:,it+1);
        J_ant = J;
        d = dmax;
        d_clk=dmax_clk;
        d_vel = dmax_vel;
        d_clkDrift=dmax_clkDrift;
    else
        gamma_est(:,it+1) = gamma_est(:,it);
        gamma_vel(:,it+1) = gamma_vel(:,it);

        ExtRxClkDrift(:,it+1)=ExtRxClkDrift(:,it);
        EstRxClkBias(:,it+1)=EstRxClkBias(:,it);
        amp_est(:,it+1) = amp_est(:,it);
        d = d/contraction;
        d_clk=d_clk/contraction;
        d_vel=d_vel/contraction;
        d_clkDrift=d_clkDrift/contraction;
    end

    if d < dmin
        d = dmax;
    end
    
    if d_clk < dmin_clk
        d_clk =dmax_clk;
    end

    if d_vel < dmin_vel
        d_vel =dmax_vel;
    end

    if d_clkDrift < dmin_clkDrift
        d_clkDrift =dmax_clkDrift;
    end

end

[navSolutions.DPE_latitude(currMeasNr), ...
 navSolutions.DPE_longitude(currMeasNr), ...
 navSolutions.DPE_height(currMeasNr)] = cart2geo(...
                                    gamma_est(1,it+1), ...
                                    gamma_est(2,it+1), ...
                                    gamma_est(3,it+1), ...
                                    5);

% Get PF-DPE estimate
DPE_X       = gamma_est(1,it+1);
DPE_Y       = gamma_est(2,it+1);
DPE_Z       = gamma_est(3,it+1);
DPE_clkBias = EstRxClkBias(:,it+1);
DPE_Xvel    = gamma_vel(1,it+1);
DPE_Yvel    = gamma_vel(2,it+1);
DPE_Zvel    = gamma_vel(3,it+1);
DPE_clkDrift= ExtRxClkDrift(:,it+1);


% === Record DPE estimates ================================================
navSolutions.DPE_estimate(currMeasNr, 1:8) =  ...
    [DPE_X DPE_Y DPE_Z DPE_clkBias DPE_Xvel DPE_Yvel DPE_Zvel DPE_clkDrift];


end


