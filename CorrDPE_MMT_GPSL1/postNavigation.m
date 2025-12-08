function [navSolutions, eph] = postNavigation(trackResults, settings, fid)
%Function calculates navigation solutions for the receiver (pseudoranges,
%positions). At the end it converts coordinates from the WGS84 system to
%the UTM, geocentric or any additional coordinate system.
%
%[navSolutions, eph] = postNavigation(trackResults, settings)
%
%   Inputs:
%       trackResults    - results from the tracking function (structure
%                       array).
%       settings        - receiver settings.
%       fid             - raw signal identifier
%   Outputs:
%       navSolutions    - contains measured pseudoranges, receiver
%                       clock error, receiver coordinates in several
%                       coordinate systems (at least ECEF and UTM).
%       eph             - received ephemerides of all SV (structure array).

%--------------------------------------------------------------------------
%                           SoftGNSS v3.0
% 
% Copyright (C) Darius Plausinaitis
% Written by Darius Plausinaitis with help from Kristin Larson
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
%$Id: postNavigation.m,v 1.1.2.22 2006/08/09 17:20:11 dpl Exp $

%% Check is there enough data to obtain any navigation solution ===========
% It is necessary to have at least three subframes (number 1, 2 and 3) to
% find satellite coordinates. Then receiver position can be found too.
% The function requires all 5 subframes, because the tracking starts at
% arbitrary point. Therefore the first received subframes can be any three
% from the 5.
% One subframe length is 6 seconds, therefore we need at least 30 sec long
% record (5 * 6 = 30 sec = 30000ms). We add extra seconds for the cases,
% when tracking has started in a middle of a subframe.

if (settings.msToProcess < 36000) 
    % Show the error message and exit
    disp('Record is to short. Exiting!');
    navSolutions = [];
    eph          = [];
    return
end

%% Pre-allocate space =====================================================
% Starting positions of the first message in the input bit stream 
% trackResults.I_P in each channel. The position is PRN code count
% since start of tracking. Corresponding value will be set to inf 
% if no valid preambles were detected in the channel.
subFrameStart  = inf(1, settings.numberOfChannels);

% Time Of Week (TOW) of the first message(in seconds). Corresponding value
% will be set to inf if no valid preambles were detected in the channel.
TOW  = inf(1, settings.numberOfChannels);

%--- Make a list of channels excluding not tracking channels --------------
activeChnList = find([trackResults.status] ~= '-');

%% Decode ephemerides =====================================================
for channelNr = activeChnList
    
    % Get PRN of current channel
    PRN = trackResults(channelNr).PRN;
    
    fprintf('Decoding NAV for PRN %02d -------------------- \n', PRN);
    %=== Decode ephemerides and TOW of the first sub-frame ================
    try
    % Edited by Sergio Vicenzo to add settings as input of NAVdecoding
    % 06 Feb 2025
    [eph(PRN), subFrameStart(channelNr), TOW(channelNr)] = ...
                                  NAVdecoding(trackResults(channelNr).I_P,...
                                  settings);  %#ok<AGROW>
    catch
    end

    try
    %--- Exclude satellite if it does not have the necessary nav data -----
    if (isempty(eph(PRN).IODC) || isempty(eph(PRN).IODE_sf2) || ...
        isempty(eph(PRN).IODE_sf3))

        %--- Exclude channel from the list (from further processing) ------
        activeChnList = setdiff(activeChnList, channelNr);
        fprintf('    Ephemeris decoding fails for PRN %02d !\n', PRN);
    else
        fprintf('    Three requisite messages for PRN %02d all decoded!\n', PRN);
    end
    catch
        %--- Exclude channel from the list (from further processing) ------
        activeChnList = setdiff(activeChnList, channelNr);
        fprintf('    Ephemeris decoding fails for PRN %02d !\n', PRN);
    end
end

%% Check if the number of satellites is still above 3 =====================
if (isempty(activeChnList) || (size(activeChnList, 2) < 4))
    % Show error message and exit
    disp('Too few satellites with ephemeris data for postion calculations. Exiting!');
    navSolutions = [];
    eph          = [];
    return
end

%% Set measurement-time point and step  ===================================
% Find start and end of measurement point locations in IF signal stream with available
% measurements
sampleStart = zeros(1, settings.numberOfChannels);
sampleEnd = inf(1, settings.numberOfChannels);

for channelNr = activeChnList
    sampleStart(channelNr) = ...
          trackResults(channelNr).absoluteSample(subFrameStart(channelNr));
    
    sampleEnd(channelNr) = trackResults(channelNr).absoluteSample(end);
end

% Second term is to make space to aviod index exceeds matrix dimensions, 
% thus a margin of 1 is added.
sampleStart = max(sampleStart) + 1;  
sampleEnd = min(sampleEnd) - 1;

%--- Measurement step in unit of IF samples -------------------------------
measSampleStep = fix(settings.samplingFreq * settings.navSolPeriod/1000);

%---  Number of measurment point from measurment start to end ------------- 
measNrSum = fix((sampleEnd-sampleStart)/measSampleStep);

%% Initialization =========================================================
% Set the satellite elevations array to INF to include all satellites for
% the first calculation of receiver position. There is no reference point
% to find the elevation angle as there is no receiver position estimate at
% this point.
satElev  = inf(1, settings.numberOfChannels);

% Save the active channel list. The list contains satellites that are
% tracked and have the required ephemeris data. In the next step the list
% will depend on each satellite's elevation angle, which will change over
% time.  
readyChnList = activeChnList;

% Set local time to inf for first calculation of receiver position. After
% first fix, localTime will be updated by measurement sample step.
localTime = inf;

%##########################################################################
%#   Do the satellite and receiver position calculations                  #
%##########################################################################

fprintf('Positions are being computed. Please wait... \n');
for currMeasNr = 1:measNrSum

    fprintf('Fix: Processing %02d of %02d \n', currMeasNr,measNrSum);
    %% Initialization of current measurement ==============================          
    % Exclude satellites, that are belove elevation mask 
    activeChnList = intersect(find(satElev >= settings.elevationMask), ...
                              readyChnList);

    % Save list of satellites used for position calculation
    navSolutions.PRN(activeChnList, currMeasNr) = ...
                                        [trackResults(activeChnList).PRN]; 

    % These two lines help the skyPlot function. The satellites excluded
    % do to elevation mask will not "jump" to possition (0,0) in the sky
    % plot.
    navSolutions.el(:, currMeasNr) = NaN(settings.numberOfChannels, 1);
    navSolutions.az(:, currMeasNr) = NaN(settings.numberOfChannels, 1);
                                     
    % Signal transmitting time of each channel at measurement sample location
    navSolutions.transmitTime(:, currMeasNr) = ...
                                         NaN(settings.numberOfChannels, 1);
    navSolutions.satClkCorr(:, currMeasNr) = ...
                                         NaN(settings.numberOfChannels, 1);                                                                  
       
    % Position index of current measurement time in IF signal stream
    % (in unit IF signal sample point)
    currMeasSample = sampleStart + measSampleStep*(currMeasNr-1);
    
%% Find pseudoranges ======================================================
    % Raw pseudorange = (localTime - transmitTime) * light speed (in m)
    % All output are 1 by settings.numberOfChannels columme vecters.
    [navSolutions.rawP(:, currMeasNr),transmitTime,localTime,codePhase]=  ...
                     calculatePseudoranges(trackResults,subFrameStart,TOW, ...
                     currMeasSample,localTime,activeChnList, settings);
    
    % Save transmitTime
    navSolutions.transmitTime(activeChnList, currMeasNr) = ...
                                        transmitTime(activeChnList);

%% Find satellites positions and clocks corrections =======================
    % Outputs are all colume vectors corresponding to activeChnList
    [satPositions, satClkCorr] = satpos(transmitTime(activeChnList), ...
                                 [trackResults(activeChnList).PRN], eph); 
                                    
     
    % Save satClkCorr
    navSolutions.satClkCorr(activeChnList, currMeasNr) = satClkCorr;

%% Find receiver position =================================================
    % 3D receiver position can be found only if signals from more than 3
    % satellites are available  
    if size(activeChnList, 2) > 3

        %=== Calculate receiver position ==================================
        % Correct pseudorange for SV clock error
        clkCorrRawP = navSolutions.rawP(activeChnList, currMeasNr)' + ...
                                                   satClkCorr * settings.c;

        % Calculate receiver position
        % Edited to also output latest satellite position for DPE
        % Sergio Vicenzo - 16 Feb 2024
        [xyzdt,navSolutions.el(activeChnList, currMeasNr), ...
               navSolutions.az(activeChnList, currMeasNr), ...
               navSolutions.DOP(:, currMeasNr),satPositions] =...
                       leastSquarePos(satPositions, clkCorrRawP, settings);

        %=== Save results =================================================
        % Receiver position in ECEF
        navSolutions.X(currMeasNr)  = xyzdt(1);
        navSolutions.Y(currMeasNr)  = xyzdt(2);
        navSolutions.Z(currMeasNr)  = xyzdt(3);       
        
		% For first calculation of solution, clock error will be set 
        % to be zero
        if (currMeasNr == 1)
        navSolutions.dt(currMeasNr) = 0;  % in unit of (m)
        else
            navSolutions.dt(currMeasNr) = xyzdt(4);  
        end
                
      
        % Save current measurement sample location 
        navSolutions.currMeasSample(currMeasNr) = currMeasSample;

        % Update the satellites elevations vector
        satElev = navSolutions.el(:, currMeasNr)';

        %=== Correct pseudorange measurements for clocks errors ===========
        navSolutions.correctedP(activeChnList, currMeasNr) = ...
                navSolutions.rawP(activeChnList, currMeasNr) + ...
                satClkCorr' * settings.c - xyzdt(4);

%% Coordinate conversion ==================================================

        %=== Convert to geodetic coordinates ==============================
        [navSolutions.latitude(currMeasNr), ...
         navSolutions.longitude(currMeasNr), ...
         navSolutions.height(currMeasNr)] = cart2geo(...
                                            navSolutions.X(currMeasNr), ...
                                            navSolutions.Y(currMeasNr), ...
                                            navSolutions.Z(currMeasNr), ...
                                            5);
        
        %=== Convert to UTM coordinate system =============================
        navSolutions.utmZone = findUtmZone(navSolutions.latitude(currMeasNr), ...
                                       navSolutions.longitude(currMeasNr));
        
        % Position in ENU
        [navSolutions.E(currMeasNr), ...
         navSolutions.N(currMeasNr), ...
         navSolutions.U(currMeasNr)] = cart2utm(xyzdt(1), xyzdt(2), ...
                                                xyzdt(3), ...
                                                navSolutions.utmZone);
        
    else
        %--- There are not enough satellites to find 3D position ----------
        disp(['   Measurement No. ', num2str(currMeasNr), ...
                       ': Not enough information for position solution.']);

        %--- Set the missing solutions to NaN. These results will be
        %excluded automatically in all plots. For DOP it is easier to use
        %zeros. NaN values might need to be excluded from results in some
        %of further processing to obtain correct results.
        navSolutions.X(currMeasNr)           = NaN;
        navSolutions.Y(currMeasNr)           = NaN;
        navSolutions.Z(currMeasNr)           = NaN;
        navSolutions.dt(currMeasNr)          = NaN;
        navSolutions.DOP(:, currMeasNr)      = zeros(5, 1);
        navSolutions.latitude(currMeasNr)    = NaN;
        navSolutions.longitude(currMeasNr)   = NaN;
        navSolutions.height(currMeasNr)      = NaN;
        navSolutions.E(currMeasNr)           = NaN;
        navSolutions.N(currMeasNr)           = NaN;
        navSolutions.U(currMeasNr)           = NaN;

        navSolutions.az(activeChnList, currMeasNr) = ...
                                             NaN(1, length(activeChnList));
        navSolutions.el(activeChnList, currMeasNr) = ...
                                             NaN(1, length(activeChnList));

        % TODO: Know issue. Satellite positions are not updated if the
        % satellites are excluded do to elevation mask. Therefore rising
        % satellites will not be included even if they will be above
        % elevation mask at some point. This would be a good place to
        % update positions of the excluded satellites.

    end % if size(activeChnList, 2) > 3

 
    %=== Direct Position Estimation plug-in module ========================
    % Added by Sergio Vicenzo - 12 Feb 2024
    tic
    [navSolutions] = DPE_module...
        (currMeasNr,navSolutions,activeChnList,...
        trackResults,currMeasSample,satPositions,...
        transmitTime(activeChnList),localTime,...
        settings,satElev,fid,xyzdt(4),satClkCorr);
    % Record DPE's processing time
    navSolutions.DPE_processingtime(currMeasNr) = toc;

    %=== Correct local time by clock error estimation =====================
    % Moved here so that localTime provided to the DPE module is
    % uncorrected by LS clock bias estimates
    % Sergio Vicenzo - 14 Feb 2024
    localTime = localTime - xyzdt(4)/settings.c;       
    navSolutions.localTime(currMeasNr) = localTime;

    %=== Update local time by measurement  step ===========================
    localTime = localTime + measSampleStep/settings.samplingFreq ;

    %=== Obtain 2D positioning error from both STL and DPE ================
    % Added by Sergio Vicenzo - 13 Feb 2024
    m2lat = 1/110734;
    m2lon = 1/103043;

    navSolutions.LLH_error(currMeasNr,1)=...
        norm(([navSolutions.latitude(currMeasNr),...
        navSolutions.longitude(currMeasNr)]...
        -settings.gt_llh(1:2))./[m2lat m2lon]);

    navSolutions.LLH_error(currMeasNr,2)=...
    norm(([navSolutions.DPE_latitude(currMeasNr),...
        navSolutions.DPE_longitude(currMeasNr)]...
        -settings.gt_llh(1:2))./[m2lat m2lon]);

%     === Prints the 2D errors of both Least Squares and DPE ==============
%     Added by Sergio Vicenzo - 14 Feb 2024
    fprintf('\nCurrent 2D Error of DPE     : %1f\n',...
        (navSolutions.LLH_error(currMeasNr,2)));
    fprintf('Current 2D Error of LS      : %1f\n',...
        (navSolutions.LLH_error(currMeasNr,1)));

    % === Calculate 3D errors =============================================
    % Added on 16 July 2024
    GT_ECEF =  llh2xyz(settings.gt_llh.*[pi/180 pi/180 1]);
    navSolutions.LLH_error(currMeasNr,3)= ...
        sqrt(sum(([navSolutions.X(currMeasNr)...
        navSolutions.Y(currMeasNr)...
        navSolutions.Z(currMeasNr)]-GT_ECEF).^2));

    pos_xyz = ...
        llh2xyz([navSolutions.DPE_latitude(currMeasNr)/180*pi,...
        navSolutions.DPE_longitude(currMeasNr)/180*pi,...
        navSolutions.DPE_height(currMeasNr)]);
    
    navSolutions.LLH_error(currMeasNr,4)=  ...
        sqrt(sum(([pos_xyz]-GT_ECEF).^2));

    % === Prints the 3D errors of both Least Squares and DPE ==============

    fprintf('\nCurrent 3D Error of DPE     : %1f\n',...
        (navSolutions.LLH_error(currMeasNr,4)));
    fprintf('Current 3D Error of LS      : %1f\n',...
        (navSolutions.LLH_error(currMeasNr,3)));

end %for currMeasNr...

% === Saves the navigation solution =======================================
% Added by Sergio Vicenzo - 09 Jan 2025
today_date= strrep(datestr(now), ':', '-');
if settings.MMT == 1
    save([settings.outfile_root, '\navSol_',settings.outfile_root,...
        '_MMT_',today_date],...
        'navSolutions','settings'); 
else
    save([settings.outfile_root, '\navSol_',settings.outfile_root,...
    '_',num2str(settings.dllCorrelatorSpacing),'_',today_date,'.mat'],...
    'navSolutions','settings'); 
end

% === Plot positioning results ============================================
% Added by Sergio Vicenzo - 15 Feb 2024
figure;
dx = geoplot(navSolutions.latitude,navSolutions.longitude);
dx.Marker = 'o';
dx.MarkerFaceColor = '#E55709';
dx.MarkerEdgeColor = '#000000';
dx.LineStyle = 'none';
dx.MarkerSize = 8;
hold on

rx = geoplot(navSolutions.DPE_latitude,navSolutions.DPE_longitude);
rx.MarkerSize = 8;
rx.LineStyle = 'none';
rx.Marker = 'o';
rx.MarkerFaceColor = '#0094FF';
rx.MarkerEdgeColor = '#000000';

geoplot(settings.gt_llh(1),settings.gt_llh(2),'k.','MarkerSize',30);
geobasemap streets;
title([settings.outfile_root, ' Geoplot'],'Interpreter','none')
legend('Least Squares (Scalar Tracking)','DPE','GT');

saveas(gcf, [pwd ,'\',settings.outfile_root,'\',...
    settings.outfile_root,'_geoplot.png']);
saveas(gcf, [pwd ,'\',settings.outfile_root,'\',...
    settings.outfile_root,'_geoplot']);

% === Plot 2D positioning error ===========================================
% Added by Sergio Vicenzo - 15 Feb 2024
figure;
ab=plot(navSolutions.LLH_error(:,1));
ab.Color = '#501D8A';
ab.LineWidth = 3;

hold on;
bc=plot(navSolutions.LLH_error(:,2));
bc.Color = '#FF9200';
bc.LineWidth = 3;
legend('Least Squares (Scalar Tracking)','DPE','FontSize',12);
label_xaxis = ['Epoch (',num2str(settings.navSolPeriod),' ms)'];
xlabel(label_xaxis,'FontSize',12);
ylabel('2D error (m)','FontSize',12);
xlim([1 length(navSolutions.LLH_error(:,2))]);
h = get(gca,'Children');
set(gca,'Children',[h(2) h(1)]);
title([settings.outfile_root,' 2D Positioning Error'],...
    'Interpreter','none');
saveas(gcf, ...
    [pwd ,'\',settings.outfile_root,'\',...
    settings.outfile_root,'_2DPositioning Error']);
saveas(gcf, ...
    [pwd ,'\',settings.outfile_root,'\',...
    settings.outfile_root,'_2DPositioning Error.png']);

% === Plot 3D positioning error ===========================================
% Added by Sergio Vicenzo - 18 Oct 2024
figure;
ab=plot(navSolutions.LLH_error(:,3));
ab.Color = '#501D8A';
ab.LineWidth = 3;

hold on;
bc=plot(navSolutions.LLH_error(:,4));
bc.Color = '#FF9200';
bc.LineWidth = 3;
legend('Least Squares (Scalar Tracking)','DPE','FontSize',12);
label_xaxis = ['Epoch (',num2str(settings.navSolPeriod),' ms)'];
xlabel(label_xaxis,'FontSize',12);
ylabel('3D error (m)','FontSize',12);
xlim([1 length(navSolutions.LLH_error(:,2))]);
h = get(gca,'Children');
set(gca,'Children',[h(2) h(1)]);
title([settings.outfile_root,' 3D Positioning Error'],...
    'Interpreter','none');
saveas(gcf, ...
    [pwd ,'\',settings.outfile_root,'\',...
    settings.outfile_root,'_3DPositioning Error']);
saveas(gcf, ...
    [pwd ,'\',settings.outfile_root,'\',...
    settings.outfile_root,'_3DPositioning Error.png']);
