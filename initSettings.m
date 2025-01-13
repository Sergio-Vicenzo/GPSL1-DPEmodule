function settings = initSettings()
%Functions initializes and saves settings. Settings can be edited inside of
%the function, updated from the command line or updated using a dedicated
%GUI - "setSettings".  
%
%All settings are described inside function code.
%
%settings = initSettings()
%
%   Inputs: none
%
%   Outputs:
%       settings     - Receiver settings (a structure). 

%--------------------------------------------------------------------------
%                           SoftGNSS v3.0
% 
% Copyright (C) Darius Plausinaitis
% Written by Darius Plausinaitis
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

% CVS record:
% $Id: initSettings.m,v 1.9.2.31 2006/08/18 11:41:57 dpl Exp $

%% Processing settings ====================================================
% Number of milliseconds to be processed used 36000 + any transients (see
% below - in Nav parameters) to ensure nav subframes are provided

% Leave around 200 to 300 milliseconds of spare data length
% i.e., deduct the actual full length by 200 or 300 ms
% This is necessary for locating bit transition in tracking.m
settings.msToProcess        = 110000;        %[ms]

% Number of channels to be used for signal processing
settings.numberOfChannels   = 12;

% Move the starting point of processing. Can be used to start the signal
% processing at any point in the data record (e.g. for long records). fseek
% function is used to move the file read point, therefore advance is byte
% based only. 
settings.skipNumberOfBytes     = 0;

%% Raw signal file name and other parameter ===============================
% This is a "default" name of the data file (signal record) to be used in
% the post-processing mode
% settings.fileName           = 'F:\B210 DATA\12Nov2022_B210_4.dat';
% settings.fileName           = 'F:\B210 DATA\12Nov2022_B210_9.dat';
% settings.fileName           = 'F:\B210 DATA\11Nov2022_B210_2.dat';
% settings.fileName           = 'F:\B210 DATA\21Dec2022_B210_5.dat';
% settings.fileName           = 'F:\B210 DATA\9Nov2022_B210_1.dat';
% settings.fileName           = 'F:\B210 DATA\30Dec2022_B210_3.dat';
% settings.fileName           = 'F:\B210 DATA\7Dec2022_B210_3.dat';
% settings.fileName           = 'F:\B210 DATA\11Nov2022_B210_3.dat';
% settings.fileName           = 'F:\B210 DATA\21Dec2022_B210_5.dat';
% settings.fileName           = 'F:\B210 DATA\29Dec2022_B210_3.dat';
% settings.fileName           = 'F:\B210 DATA\29Dec2022_B210_11.dat';
settings.fileName           = 'F:\NSL STEREO\D_20191230_TSTEb_GPS_1713.dat';
% settings.fileName           = 'F:\NSL STEREO\D_20190414_YouMaTai_1920.dat';
% settings.fileName           = 'F:\NSL STEREO\d_0850_3.dat';
% settings.fileName           = 'E:\test3.dat';
% settings.fileName           ='F:\Spirent simulation data\003\File_003.bin';%'C:\Users\gnss\Desktop\Discussion\Amungo Navigation for NUT4NT+\dump_ch1.bin';  
% settings.fileName           = 'F:\N210 DATA\1Mar2024_N210_GPSL1_1.dat';
% settings.fileName           = 'F:\LABSAT\HaiDi_File_039.bin';
% Data type used to store one sample
settings.dataType           = 'schar';  % "schar" = int8, "short" = int16
 
% File Types
%1 - 8 bit real samples S0,S1,S2,...
%2 - 8 bit I/Q samples I0,Q0,I1,Q1,I2,Q2,...                      
settings.fileType           = 2;

% Intermediate, sampling and code frequencies
settings.IF                 = 0;     % [Hz]
settings.samplingFreq       = 26e6;%26e6;%58e6;        % [Hz]
settings.codeFreqBasis      = 1.023e6;     % [Hz]

% Define number of chips in a code period
settings.codeLength         = 1023;

%% Acquisition settings ===================================================
% Skips acquisition in the script postProcessing.m if set to 1
settings.skipAcquisition    = 0;
% List of satellites to look for. Some satellites can be excluded to speed
% up acquisition
settings.acqSatelliteList   = 1:32;         %[PRN numbers]
% Band around IF to search for satellite signal. Depends on max Doppler.
% It is single sideband, so the whole search band is tiwce of it.
settings.acqSearchBand      = 7000;           %[Hz]
% Threshold for the signal presence decision rule
settings.acqThreshold       = 1.8; % Original: 1.8
% Sampling rate threshold for downsampling 
settings.resamplingThreshold    = 8e6;            % [Hz]
% Enable/dissable use of downsampling for acquisition
settings.resamplingflag         = 0;              % 0 - Off
                                                  % 1 - On
%% Tracking loops settings ================================================
% Code tracking loop parameters
settings.dllDampingRatio         = 0.707; 
% settings.dllNoiseBandwidth       = 1.5;       %[Hz] 
% Disabled by Sergio Vicenzo - 06 Jan 2025
% settings.dllNoiseBandwidth is fixed at 1.5 for PDI = 1 ms 
% and 4 for PDI >= 10 ms
settings.dllCorrelatorSpacing    = 0.3;     %[chips] % prev. 0.5 

% Carrier tracking loop parameters
settings.pllDampingRatio         = 0.707; 
% settings.pllNoiseBandwidth       = 20;      %[Hz] 
% Disabled by Sergio Vicenzo - 06 Jan 2025
% settings.pllNoiseBandwidth is fixed at 20 for PDI = 1 ms 
% and 5 for PDI >= 5 ms

% Integration time for DLL and PLL
% settings.intTime                 = 0.001;      %[s] 
% Disabled by Sergio Vicenzo - 06 Jan 2025
% To change tracking coherent integration time, use settings.DPE_cohInt
% below

%% Navigation solution settings ===========================================

% Period for calculating pseudoranges and position
settings.navSolPeriod       = 500;          %[ms]

% Elevation mask to exclude signals from satellites at low elevation
settings.elevationMask      = 0;           %[degrees 0 - 90]
% Enable/dissable use of tropospheric correction
settings.useTropCorr        = 1;            % 0 - Off
                                            % 1 - On

% True position of the antenna in UTM system (if known). Otherwise enter
% all NaN's and mean position will be used as a reference .
settings.truePosition.E     = nan;
settings.truePosition.N     = nan;
settings.truePosition.U     = nan;

%% Plot settings ==========================================================
% Enable/disable plotting of the tracking results for each channel
settings.plotTracking       = 0;            % 0 - Off
                                            % 1 - On
%% Constants ==============================================================
settings.c                  = 299792458;    % The speed of light, [m/s]
settings.startOffset        = 68.802;       %[ms] Initial sign. travel time

%% CNo Settings============================================================
% Accumulation interval in Tracking (in Sec)
settings.CNo.accTime=0.001;
% Accumulation interval for computing VSM C/No (in ms)
settings.CNo.VSMinterval = 40;

%% Ground Truth ===========================================================
% For evaluation of SDR performance
% Added by Sergio Vicenzo - 13 Feb 2024
% settings.gt_llh = [31.26529162 121.62553832  0]; % File_002
% settings.gt_llh = [31.264469 121.625567 0]; % File_003
% settings.gt_llh = [22.30465290 114.18096270 10.667]; %12Nov2022_B210_4
% settings.gt_llh =[22.30396900 114.18027170 10.652]; %09 NOV 2022 (FILE 1 - 5)
% settings.gt_llh =[22.29949700 114.17810610 2.813]; %11 NOV 2022
% settings.gt_llh = [22.3122139, 114.17280277777777, 65]; %D_20190414_YouMaTai_1920
settings.gt_llh = [22.299866 114.180036 3]; % D_20191230_TSTEb_GPS_1713
% settings.gt_llh = [22.30477703, 114.18030159, 11.591]; %30 DEC 2022
% settings.gt_llh = [22.30290170 114.17775590 11.611]; %21 DEC 2022	(FILE 1 - 5)
% settings.gt_llh =[22.30729456 114.17950917 28.894]; %09 NOV 2022 (FILE 6 - 10)
% settings.gt_llh = [22.30459780 114.18011900 61.095]; %- R core rooftop
% settings.gt_llh = [22.3028020248354, 114.193081317790, 7.487]; %Whampoa
% settings.gt_llh = [22.30401360 114.17890410 9.866]; %- 12Nov2022_B210_9
% settings.gt_llh =  [22.30290170 114.17775590 9.769]; % 21 DEC 2022	(FILE 1 - 5)
% settings.gt_llh = [22.30386634 114.18001750 11.611]; %29 DEC 2022 (FILE 1 - 3)
% settings.gt_llh = [22.30536039 114.18041400 11.611]; %29 DEC 2022 (FILE 4 - 10)
% settings.gt_llh = [22.30554810 114.18013600 20.664] ; %29 DEC 2022 (FILE 11 - 15)
% settings.gt_llh = [22.30401360 114.17890410 9.866]; %7Dec2022_B210_3
% settings.gt_llh = [22.29949700 114.17810610 2.813] ; %11Nov2022_B210_2
% settings.gt_llh = [22.292958333333335 114.1742888888889 2]; %HaiDi files

%% DPE config =============================================================
% Chip spacing for DPE "pre-calculate correlations" implementation

% A traditional/Pure DPE implementation is computationally extensive. A
% "pre-calculate correlations" implementation calculates the correlations 
% per every "settings.chipspacing_dpe_precalc" number of chip spacing. 

% The correlations for each candidate position would then be given based
% on its corresponding code phase. Interpolation is used to obtain the 
% correlations for every candidate position.

% The lesser the chip spacing, the more accurate the correlations for
% each candidate position, but comes with higher computational load.

% Higher chip spacing has lower computational load, but correlations for
% each candidate position would be less accurate and rely more on 
% interpolation 

% Check "DPE_module.m" for more details

% Added by Sergio Vicenzo - 14 Feb 2024
settings.chipspacing_dpe_precalc = settings.codeFreqBasis ...
                                    / settings.samplingFreq;

% DPE's non-coherent integration time
% Added by Sergio Vicenzo - 5 Mar 2024
settings.DPE_nonCohInt = 1; % in ms 
% Note that this does not change the tracking's non-coherent integration time
% Tracking non-coherent integration time is fixed at 1 

% DPE's coherent integration time
% Added by Sergio Vicenzo - 06 Jan 2025
settings.DPE_cohInt = 20; % in ms % For both STL and DPE. 
% Recommended to use 20 ms of coherent integration

% ONLY 1 ms, 10 ms, and 20 ms coherent integration AVAILABLE! 
% USING OTHER PDIs requires further modification to NAVdecoding.m

% If PDI = 10 ms fails to decode the navigation message, 
% change searchStartOffset in NAVdecoding.m

% Span of lat-long search space i.e., plus-minus
% "settings.DPE_latlong_span" meters
% Added by Sergio Vicenzo - 16 Feb 2024
settings.DPE_latlong_span = 30;

% Span of height search space i.e., plus-minus "settings.DPE_height_span"
% meters
% Large search space for height estimation due to GNSS' inherent 
% nature of high uncertainty in height estimation
% Added by Sergio Vicenzo - 16 Feb 2024
settings.DPE_height_span = 50;

% Span of clock bias search space i.e., plus-minus 
% "settings.DPE_clkBias_span" meters
% Added by Sergio Vicenzo - 16 Feb 2024
settings.DPE_clkBias_span = 20;

% The spacing between each candidate lat-long-height 
% "settings.candPVT_spacing" meters between the generated grid of candidate
% latitude-longitude-height
% Added by Sergio Vicenzo - 18 April 2024
settings.candPVT_spacing = 1; % in meters

% Enable/disable plotting correlogram of DPE at every epoch
% Added by Sergio Vicenzo - 5 Mar 2024
settings.DPE_plotCorrelogram = 0;   % 0 - Off
                                    % 1 - On

%% MMT config =============================================================
% Multipath Mitigation Technology to improve performance of DPE in urban
% environments
settings.MMT = 0;                   % 0 - Off
                                    % 1 - On

% Amplitudes of secondary path i.e., reflected or non-line-of-sight path
% are typically smaller than that of the line-og-sight/direct path. MMT
% performance can be improved by constraining the amplitude of NLOS path to
% be lower than the LOS by settings.MMT_const times

% If settings.MMT_const  = 1, reflected path amplitude is assumed to be 
% able to have the same amplitude as the LOS path

% Added by Sergio Vicenzo - 15 Mar 2024
settings.MMT_const = 0.8; 
