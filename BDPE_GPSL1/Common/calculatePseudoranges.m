function [pseudoranges,transmitTime,localTime,codePhase] = ...
             calculatePseudoranges(trackResults,subFrameStart,TOW,currMeasSample, ...
             localTime,channelList, settings)
         
%calculatePseudoranges finds relative pseudoranges for all satellites
%listed in CHANNELLIST at the specified millisecond of the processed
%signal. The pseudoranges contain unknown receiver clock offset. It can be
%found by the least squares position search procedure. 
%
% Edited to also output code phase measurements for DPE plug-in module
% Sergio Vicenzo - 14 Feb 2024

% [pseudoranges,transmitTime,localTime] = ...
%              calculatePseudoranges(trackResults,subFrameStart,TOW,currMeasSample, ...
%              localTime,channelList, settings)
%
%   Inputs:
%       trackResults    - output from the tracking function
%       subFrameStart   - the array contains positions of the first
%                       preamble in each channel. The position is ms count 
%                       since start of tracking. Corresponding value will
%                       be set to 0 if no valid preambles were detected in
%                       the channel: 
%                       1 by settings.numberOfChannels
%       TOW             - Time Of Week (TOW) of the first sub-frame in the bit
%                       stream (in seconds)
%       currMeasSample  - current measurement sample location(measurement time)
%       localTime       - local time(in GPST) at measurement time
%       channelList     - list of channels to be processed
%       settings        - receiver settings
%
%   Outputs:
%       pseudoranges    - relative pseudoranges to the satellites. 

%       transmitTime    - transmitting time of channels to be processed 
%                         corresponding to measurement time  
%       localTime       - local time(in GPST) at measurement time
%       codePhase       - Code phase from start of a PRN code 
%                         to current measurement sample location for DPE
%                       

%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%                           SoftGNSS v3.0
% 
% Copyright (C) Darius Plausinaitis
% Written by Yafeng Li
% Based on Darius Plausinaitis, Peter Rinder and Nicolaj Bertelsen
%--------------------------------------------------------------------------
%
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
% $Id: calculatePseudoranges.m,v 1.1.2.18 2006/08/09 17:20:11 dpl Exp $
    
% Transmitting Time of all channels at current measurement sample location
transmitTime = inf(1, settings.numberOfChannels);

% Initialize parameters
% Added by Sergio Vicenzo - 19 Feb 2024
codePhase = zeros(length(channelList),1);

%--- For all channels in the list ... 
for channelNr = channelList
    
    % Find index of I_P stream whose integration contains current 
    % measurment point location 
    for index = 1: length(trackResults(channelNr).absoluteSample)
        if(trackResults(channelNr).absoluteSample(index) > currMeasSample )
            break
        end 
    end
    index = index - 1;

          
    % Update the phasestep based on code freq and sampling frequency
    codePhaseStep = trackResults(channelNr).codeFreq(index) / settings.samplingFreq ;
    
    % Code phase from start of a PRN code to current measurement sample location 
    % Edited to record code phase for each satellites 
    % Sergio Vicenzo - 12 Feb 2024
    codePhase(channelNr)  = (trackResults(channelNr).remCodePhase(index) +  ...
                          codePhaseStep * (currMeasSample - ...
                          trackResults(channelNr).absoluteSample(index) ));

    
    % Transmitting Time (in unite of s)at current measurement sample location
    % codePhase/settings.codeLength: fraction part of a PRN code
    % index - subFrameStart(channelNr): integer number of PRN code
    % Edited by Sergio Vicenzo - 06 Jan 2025
    transmitTime(channelNr) =  ...
        (codePhase(channelNr)/settings.codeLength/settings.DPE_cohInt ...
                          + index  - ...
                          subFrameStart(channelNr)) * ...
                          settings.codeLength * settings.DPE_cohInt /...
                          settings.codeFreqBasis + TOW(channelNr);

end

% At first time of fix, local time is initialized by transmitTime and 
% settings.startOffset
if (localTime == inf)
    maxTime   = max(transmitTime(channelList));
    localTime = maxTime + settings.startOffset/1000;  
end


%--- Convert travel time to a distance ------------------------------------
% The speed of light must be converted from meters per second to meters
% per millisecond. 
pseudoranges    = (localTime - transmitTime) * settings.c; 
