
addpath include               % The software receiver functions
addpath Common       % Common functions between differnt SDR receivers
% addpath ('../IF_Data_Set')    % IF data sets for each SDR receivers
activeChnList = find([trackResults.status] ~= '-');

%% Decode ephemerides =====================================================
for channelNr = activeChnList
    
    % Get PRN of current channel
    PRN = trackResults(channelNr).PRN;
    
    fprintf('Decoding NAV for PRN %02d -------------------- \n', PRN);
    %=== Decode ephemerides and TOW of the first sub-frame ================
%     try
    % Edited by Sergio Vicenzo to add settings as input of NAVdecoding
    % 06 Feb 2025
    [eph(PRN), subFrameStart(channelNr), TOW(channelNr)] = ...
                                  NAVdecoding(trackResults(channelNr).I_P,...
                                  settings);  %#ok<AGROW>
%     catch
%     end

end