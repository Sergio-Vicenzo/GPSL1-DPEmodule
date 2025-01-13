function [eph, subFrameStart,TOW] = NAVdecoding(I_P_InputBits,settings)

% findPreambles finds the first preamble occurrence in the bit stream of
% each channel. The preamble is verified by check of the spacing between
% preambles (6sec) and parity checking of the first two words in a
% subframe. At the same time function returns list of channels, that are in
% tracking state and with valid preambles in the nav data stream.
%
%[eph, subFrameStart,TOW] = CNAVdecoding(I_P_InputBits)
%
%   Inputs:
%       I_P_InputBits   - output from the tracking function
%
%   Outputs:
%       subFrameStart   - Starting positions of the first message in the 
%                       input bit stream I_P_InputBits in each channel. 
%                       The position is CNAV bit(20ms before convolutional decoding) 
%                       count since start of tracking. Corresponding value will
%                       be set to inf if no valid preambles were detected in
%                       the channel.
%       TOW             - Time Of Week (TOW) of the first message(in seconds).
%                       Corresponding value will be set to inf if no valid preambles
%                       were detected in the channel.
%       eph             - SV ephemeris. 

%--------------------------------------------------------------------------
%                           SoftGNSS v3.0
% 
% Copyright (C) Yafeng Li, Denish Akos
% Written by Yafeng Li
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
% $Id: findPreambles.m,v 1.1.2.10 2017/01/19 21:13:22 dpl Exp $


%--- Initialize ephemeris structute  --------------------------------------
% This is in order to make sure variable 'eph' for each SV has a similar 
% structure when only one or even none of the three requisite messages
% is decoded for a given PRN.
eph = eph_structure_init();

% Starting positions of the first message in the input bit stream
subFrameStart = inf;

% TOW of the first message
TOW = inf;

%% Bit and frame synchronization ====================================
% Preamble search can be delayed to a later point in the tracking results
% to avoid noise due to tracking loop transients
searchStartOffset = 1000; 
% Try to change this parameter if the navigation message decoding fails
% Sergio Vicenzo - 07 Jan 2025

%--- Generate the preamble pattern ----------------------------------------
preamble_bits = [1 -1 -1 -1 1 -1 1 1];

% "Upsample" the preamble - make 20 vales per one bit. The preamble must be
% found with precision of a sample.

% Edited by Sergio Vicenzo to account various coherent integration
% tracking
% 06 Jan 2025
preamble_ms = kron(preamble_bits, ones(1, 20/settings.DPE_cohInt));

% Correlate tracking output with preamble =================================
% Read output from tracking. It contains the navigation bits. The start
% of record is skiped here to avoid tracking loop transients.
bits = I_P_InputBits(1 + searchStartOffset : end);

% Now threshold the output and convert it to -1 and +1
bits(bits > 0)  =  1;
bits(bits <= 0) = -1;

% Correlate tracking output with the preamble
tlmXcorrResult = xcorr(bits, preamble_ms);

% Find all starting points off all preamble like patterns =================
clear index
clear index2

xcorrLength = (length(tlmXcorrResult) +  1) /2;

%--- Find at what index/ms the preambles start ------------------------
% Edited by Sergio Vicenzo to account various coherent integration
% tracking
% 06 Jan 2025
if settings.DPE_cohInt == 1

    index = find(...
        abs(tlmXcorrResult(xcorrLength : xcorrLength * 2 - 1)) > 153)' + ...
        searchStartOffset;

else % Coherent integration is higher than 1 ms

    index = find(...
    abs(tlmXcorrResult(xcorrLength : xcorrLength * 2 - 1)) > 6)' + ...
    searchStartOffset;

end

% Analyze detected preamble like patterns ================================
 for i = 1:size(index) % For each occurrence
    
    %--- Find distances in time between this occurrence and the rest of
    %preambles like patterns. If the distance is 6000 milliseconds (one
    %subframe), the do further verifications by validating the parities
    %of two GPS words
    
    index2 = index - index(i);
    
    if (~isempty(find(index2 == 6000/settings.DPE_cohInt, 1)))
        
        %=== Re-read bit values for preamble verification ==============
        % Preamble occurrence is verified by checking the parity of
        % the first two words in the subframe. Now it is assumed that
        % bit boundaries a known. Therefore the bit values over 20ms are
        % combined to increase receiver performance for noisy signals.
        % in Total 62 bits mast be read :
        % 2 bits from previous subframe are needed for parity checking;
        % 60 bits for the first two 30bit words (TLM and HOW words).
        % The index is pointing at the start of TLM word.

        % Edited by Sergio Vicenzo to account various coherent integration
        % tracking
        % 06 Jan 2025
        bits = I_P_InputBits(index(i)-40/settings.DPE_cohInt...
            : (index(i) + (20 * 60) / settings.DPE_cohInt)-1)';
        
        %--- Combine the 20 values of each bit ------------------------
        % Edited by Sergio Vicenzo - 06 Jan 2025
        if settings.DPE_cohInt < 20
            bits = reshape(bits, 20 / settings.DPE_cohInt, ...
                ((size(bits, 1) / 20)* settings.DPE_cohInt));
            bits = sum(bits,1); 
        else % Coherent integration is 20 ms
            bits=bits';
        end
        
        % Now threshold and make it -1 and +1
        bits(bits > 0)  = 1;
        bits(bits <= 0) = -1;
        
        %--- Check the parity of the TLM and HOW words ----------------
        if (navPartyChk(bits(1:32)) ~= 0) && ...
                (navPartyChk(bits(31:62)) ~= 0)
            % Parity was OK. Record the preamble start position. Skip
            % the rest of preamble pattern checking for this channel
            % and process next channel.
            
            subFrameStart = index(i);
            break;
        end % if parity is OK ...
        
    end % if (~isempty(find(index2 == 6000)))
end % for i = 1:size(index)

% Exclude channel from the active channel list if no valid preamble was
% detected
if subFrameStart == inf
    disp('Could not find valid preambles in channel! ');
    return
end
    
%% Decode ephemerides ===============================================
%=== Convert tracking output to navigation bits =======================
%--- Copy 5 sub-frames long record from tracking output ---------------
% Edited by Sergio Vicenzo - 06 Jan 2025
navBitsSamples = I_P_InputBits(subFrameStart - 20/settings.DPE_cohInt: ...
    subFrameStart + (1500 * 20)/settings.DPE_cohInt -1)';

%--- Group every 20 values of bits into columns ------------------------
% Edited by Sergio Vicenzo - 06 Jan 2025
if settings.DPE_cohInt<20
    navBitsSamples = reshape(navBitsSamples, ...
        20/settings.DPE_cohInt, ...
        (size(navBitsSamples, 1) / 20)* settings.DPE_cohInt);
else
    % If PDI = 20 ms, the bits do not need to be grouped again
    navBitsSamples=navBitsSamples';
end

%--- Sum all samples in the bits to get the best estimate -------------
% Edited by Sergio Vicenzo - 06 Jan 2025
navBits = sum(navBitsSamples,1); % Does nothing when PDI = 20 ms

%--- Now threshold and make 1 and 0 -----------------------------------
% The expression (navBits > 0) returns an array with elements set to 1
% if the condition is met and set to 0 if it is not met.
navBits = (navBits > 0);

%--- Convert from decimal to binary -----------------------------------
% The function ephemeris expects input in binary form. In Matlab it is
% a string array containing only "0" and "1" characters.
navBitsBin = dec2bin(navBits);

%=== Decode ephemerides and TOW of the first sub-frame ================
[eph, TOW] = ephemeris(navBitsBin(2:1501)', navBitsBin(1));
